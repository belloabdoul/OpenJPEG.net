#region License
/*
 * Copyright (c) 2002-2007, Communications and Remote Sensing Laboratory, Universite catholique de Louvain (UCL), Belgium
 * Copyright (c) 2002-2007, Professor Benoit Macq
 * Copyright (c) 2001-2003, David Janssens
 * Copyright (c) 2002-2003, Yannick Verschueren
 * Copyright (c) 2003-2007, Francois-Olivier Devaux and Antonin Descampe
 * Copyright (c) 2005, Herve Drolon, FreeImage Team
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS `AS IS'
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */
#endregion
//STANDARD_SLOW_VERSION Uses the OpenJpeg 2.1 codepath for DWT decoding.
//#define STANDARD_SLOW_VERSION
using System;
using System.Diagnostics;
using System.Threading;

namespace OpenJpeg.Internal;

/// <summary>
/// Discrete Wavelet Transform
/// </summary>
/// <remarks>
/// The DWT can be reversible, or irreversible.
/// 
/// Reversible uses integer math and a modified YUV colorspace,
/// irreversible uses floating point math and the YCbCr colorspace
/// 
/// Irreversible is only applicable for lossy encoding.
/// </remarks>
internal static class Dwt
{

    #region Consts

    private const float K = 1.230174105f;
    private const float InvK = (float)(1.0 / 1.230174105);
    private const float C13318 = 1.625732422f; //<-- two_invK
    private const float DwtAlpha = -1.586134342f; //  12994
    private const float DwtBeta = -0.052980118f; //    434
    private const float DwtGamma = 0.882911075f; //  -7233
    private const float DwtDelta = 0.443506852f; //  -3633

    /// <summary>
    /// Number of int32 values in a SSE2 register
    /// </summary>
    /// <remarks>
    /// We don't currently support SSE2, but maybe in the future
    /// </remarks>
    private const uint VregIntCount = 4;

    private const uint ParallelCols53 = 2 * VregIntCount;
    private const int NbEltsV8 = 8;

    /// <summary>
    /// This table contains the norms of the 9-7 wavelets for different bands.
    /// </summary>
    private static readonly double[][] DwtNormsReal = {
        new double[] {1.000, 1.965, 4.177, 8.403, 16.90, 33.84, 67.69, 135.3, 270.6, 540.9},
        new double[] {2.022, 3.989, 8.355, 17.04, 34.27, 68.63, 137.3, 274.6, 549.0},
        new double[] {2.022, 3.989, 8.355, 17.04, 34.27, 68.63, 137.3, 274.6, 549.0},
        new double[] {2.080, 3.865, 8.307, 17.18, 34.71, 69.59, 139.3, 278.6, 557.2}
    };

    /// <summary>
    /// This table contains the norms of the 5-3 wavelets for different bands.
    /// </summary>
    private static readonly double[][] DwtNorms = {
        new double[] {1.000, 1.500, 2.750, 5.375, 10.68, 21.34, 42.67, 85.33, 170.7, 341.3},
        new double[] {1.038, 1.592, 2.919, 5.703, 11.33, 22.64, 45.25, 90.48, 180.9},
        new double[] {1.038, 1.592, 2.919, 5.703, 11.33, 22.64, 45.25, 90.48, 180.9},
        new double[] {.7186, .9218, 1.586, 3.043, 6.019, 12.01, 24.00, 47.97, 95.93}
    };

    #endregion

    private delegate void Encodefunc(int[] a, int dn, int sn, int cas);

    private delegate void EncodeAndDeinterleaveVfunc(int[] a, int aPt, int[] tmp, uint height, bool even, uint strideWidth, uint cols);

    private delegate void EncodeAndDeinterleaveHOneRowfunc(int[] row, int rowPt, int[] tmp, uint width, bool even);

    /// <summary>
    /// Forward 5-3 wavelet transform in 2-D
    /// </summary>
    /// <remarks>2.5 - opj_dwt_encode</remarks>
    internal static bool Encode(TileCoder tcd, TcdTilecomp tilec)
    {
        return EncodeProcedure(tilec, EncodeAndDeinterleaveV, EncodeAndDeinterleaveH_OneRow, tcd.DisableMultiThreading);
    }

    /// <summary>
    /// Forward 9-7 wavelet transform in 2-D
    /// </summary>
    /// <remarks>2.5 - opj_dwt_encode_real</remarks>
    internal static bool EncodeReal(TileCoder tcd, TcdTilecomp tilec)
    {
        return EncodeProcedure(tilec, EncodeAndDeinterleaveV_Real, EncodeAndDeinterleaveH_OneRowReal, tcd.DisableMultiThreading);
    }

    //2.5 - opj_dwt_encode_procedure
    private static bool EncodeProcedure(TcdTilecomp tilec, 
        EncodeAndDeinterleaveVfunc encodeAndDeinterleaveV, 
        EncodeAndDeinterleaveHOneRowfunc encodeAndDeinterleaveHOneRow, 
        bool disableMultiThreading)
    {
        int numThreads;
        ThreadPool.GetAvailableThreads(out numThreads, out _);
        numThreads = disableMultiThreading ? 1 : Math.Min(Environment.ProcessorCount, numThreads);

        var tiledp = tilec.data;

        var w = (uint)(tilec.x1 - tilec.x0);
        var l = (int)tilec.numresolutions - 1;

        var tr = tilec.resolutions;
        var curRes = l; //<-- pointer to tilec.resolutions
        var lastRes = curRes - 1; //<-- pointer to tilec.resolutions

        var dataSize = MaxResolution(tilec.resolutions, (int)tilec.numresolutions);
        if (dataSize > Constants.SizeMax / (NbEltsV8 * sizeof(int)))
            return false;
        dataSize *= NbEltsV8; //C# org impl is number of bytes, here it's number of ints
        var bj = new int[dataSize];
        var i = l;

        using (var reset = new ManualResetEvent(false))
        {
            while (i-- != 0)
            {
                //Width of the resolution level computed
                var rw = (uint)(tr[curRes].x1 - tr[curRes].x0);

                //Height of the resolution level computed
                var rh = (uint)(tr[curRes].y1 - tr[curRes].y0);

                //Width of the resolution level once lower than computed one
                var rw1 = (uint)(tr[lastRes].x1 - tr[lastRes].x0);

                //Height of the resolution level once lower than computed one 
                var rh1 = (uint)(tr[lastRes].y1 - tr[lastRes].y0);

                //0 = non inversion on vertical filtering 1 = inversion between low-pass and high-pass filtering
                var casRow = tr[curRes].x0 & 1;

                //0 = non inversion on horizontal filtering 1 = inversion between low-pass and high-pass filtering
                var casCol = tr[curRes].y0 & 1;

                var sn = (int)rh1;
                var dn = (int)(rh - rh1);

                // Perform vertical pass
                if (numThreads <= 1 || rw < 2 * NbEltsV8)
                {
                    var j = 0;
                    for (; j + NbEltsV8 - 1 < rw; j += NbEltsV8)
                    {
                        encodeAndDeinterleaveV(tiledp, j, bj, rh, casCol == 0, w, NbEltsV8);
                    }
                    if (j < rw)
                    {
                        encodeAndDeinterleaveV(tiledp, j, bj, rh, casCol == 0, w, rw - (uint)j);
                    }
                }
                else
                {
                    var numJobs = numThreads;

                    if (rw < numJobs)
                    {
                        numJobs = (int)rw;
                    }

                    var stepJ = rw / (uint)numJobs / NbEltsV8 * NbEltsV8;

                    reset.Reset();
                    //Alternativly, we can set this to num_jobs and remove the Interlocked.Increment
                    //and the Interlocked.Decrement after the for loop
                    var nThreadWorkers = 1;

                    for (uint j = 0; j < numJobs; j++)
                    {
                        var job = new EncodeVJob(
                            new DwtLocal() { 
                                Mem = new int[dataSize],
                                Dn = dn,
                                Sn = sn,
                                Cas = casCol
                            },
                            rh,
                            w,
                            tiledp,
                            0,
                            j * stepJ,
                            j + 1 == numJobs ? rw : (j + 1) * stepJ,
                            encodeAndDeinterleaveV
                        );

                        Interlocked.Increment(ref nThreadWorkers);
                        ThreadPool.QueueUserWorkItem((x) =>
                        {
                            try { encode_v_func((EncodeVJob)x); }
                            finally
                            {
                                if (Interlocked.Decrement(ref nThreadWorkers) == 0)
                                    reset.Set();
                            }
                        }, job);
                    }
                    if (Interlocked.Decrement(ref nThreadWorkers) == 0)
                        reset.Set();
                    reset.WaitOne();
                }

                sn = (int)rw1;
                dn = (int)(rw - rw1);

                // Perform horizontal pass
                if (numThreads <= 1 || rh <= 1)
                {
                    for (var j = 0; j < rh; j++)
                    {
                        encodeAndDeinterleaveHOneRow(tiledp, j * (int)w, bj, rw, casRow == 0);
                    }
                }
                else
                {
                    var numJobs = numThreads;

                    if (rh < numJobs)
                    {
                        numJobs = (int)rh;
                    }

                    var stepJ = rh / (uint)numJobs;

                    reset.Reset();
                    //Alternativly, we can set this to num_jobs and remove the Interlocked.Increment
                    //and the Interlocked.Decrement after the for loop
                    var nThreadWorkers = 1;

                    for (uint j = 0; j < numJobs; j++)
                    {
                        var maxJ = (j + 1) * stepJ;
                        if (j == numJobs - 1)
                            maxJ = rh;

                        var job = new EncodeHJob(
                            new DwtLocal()
                            {
                                Mem = new int[dataSize],
                                Dn = dn,
                                Sn = sn,
                                Cas = casRow
                            },
                            rw,
                            w,
                            tiledp,
                            0,
                            j * stepJ,
                            maxJ,
                            encodeAndDeinterleaveHOneRow
                        );

                        Interlocked.Increment(ref nThreadWorkers);
                        ThreadPool.QueueUserWorkItem((x) =>
                        {
                            try { encode_h_func((EncodeHJob)x); }
                            finally
                            {
                                if (Interlocked.Decrement(ref nThreadWorkers) == 0)
                                    reset.Set();
                            }
                        }, job);
                    }
                    if (Interlocked.Decrement(ref nThreadWorkers) == 0)
                        reset.Set();
                    reset.WaitOne();
                }

                curRes = lastRes;

                lastRes--;
            }
        }

        return true; ;
    }

    //2.5 - opj_dwt_encode_h_func
    private static void encode_h_func(EncodeHJob job)
    {
        for (var j = job.MinJ; j < job.MaxJ; j++)
        {
            job.Fn(job.Tiled, job.Tiledp + (int)j * (int)job.W, job.H.Mem, job.Rw, job.H.Cas == 0);
        }
    }

    //2.5 - opj_dwt_encode_v_func
    private static void encode_v_func(EncodeVJob job)
    {
        uint j;
        for (j = job.MinJ; j + NbEltsV8 - 1 < job.MaxJ; j += NbEltsV8)
        {
            job.EncodeAndDeinterleaveV(job.Tiled, job.Tiledp + (int)j, job.V.Mem, job.Rh, job.V.Cas == 0, job.W, NbEltsV8);
        }
        if (j < job.MaxJ)
        {
            job.EncodeAndDeinterleaveV(job.Tiled, job.Tiledp + (int)j, job.V.Mem, job.Rh, job.V.Cas == 0, job.W, job.MaxJ - j);
        }
    }

    /// <summary>
    /// Determine maximum computed resolution level for inverse wavelet transform
    /// </summary>
    /// <param name="r">Resolutions</param>
    /// <param name="i">Number of resolutions that will be decoded</param>
    /// <remarks>2.5 - opj_dwt_max_resolution</remarks>
    private static uint MaxResolution(TcdResolution[] rs, int numres)
    {
        uint mr = 0;
        uint w;
        for (var c = 1; c < numres; c++)
        {
            var r = rs[c];
            if (mr < (w = (uint)(r.x1 - r.x0)))
                mr = w;
            if (mr < (w = (uint)(r.y1 - r.y0)))
                mr = w;
        }
        return mr;
    }

    //2.5 - opj_dwt_calc_explicit_stepsizes
    internal static void CalcExplicitStepsizes(TileCompParams tccp, uint prec)
    {
        var numbands = 3 * tccp.numresolutions - 2;
        for (uint bandno = 0; bandno < numbands; bandno++)
        {
            double stepsize;

            var resno = bandno == 0 ? 0 : (bandno - 1) / 3 + 1;
            var orient = bandno == 0 ? 0 : (bandno - 1) % 3 + 1;
            var level = tccp.numresolutions - 1 - resno;
            var gain = tccp.qmfbid == 0 ? 0U : orient == 0 ? 0U : orient == 1 || 
                                                                  orient == 2 ? 1U : 2U;
            if (tccp.qntsty == CCP_QNTSTY.NOQNT)
                stepsize = 1.0;
            else
            {
                var norm = GetNormReal(level, orient);
                stepsize = (1 << (int)gain) / norm;
            }
            EncodeStepsize((int)Math.Floor(stepsize * 8192.0), (int)(prec + gain), out tccp.stepsizes[bandno]);
        }
    }

    //2.5 - opj_dwt_getnorm_real
    private static double GetNormReal(uint level, uint orient)
    {
        //This is just a band-aid to avoid a buffer overflow
        if (orient == 0 && level >= 10)
            level = 9;
        else if (orient > 0 && level >= 9)
            level = 8;
        return DwtNormsReal[orient][level];
    }

    /// <remarks>2.5 - opj_dwt_encode_stepsize</remarks>
    private static void EncodeStepsize(int stepsize, int numbps, out StepSize bandnoStepsize)
    {
        var p = MyMath.int_floorlog2(stepsize) - 13;
        var n = 11 - MyMath.int_floorlog2(stepsize);
        bandnoStepsize = new StepSize(
            numbps - p, 
            (n < 0 ? stepsize >> -n : stepsize << n) & 0x7ff);
    }

    //2.5 - opj_dwt_decode_real
    internal static bool DecodeReal(TileCoder tcd, TcdTilecomp tilec, uint numres)
    {
        if (tcd.WholeTileDecoding)
            return decode_tile_97(tilec, numres, tcd.DisableMultiThreading);
        else
            return decode_partial_97(tilec, numres);
    }

    //2.5 - opj_dwt_decode_partial_97
    private static bool decode_partial_97(TcdTilecomp tilec, uint numres)
    {
        var h = new V4dwtLocal();
        var v = new V4dwtLocal();
        // This value matches the maximum left/right extension given in tables
        // F.2 and F.3 of the standard. Note: in opj_tcd_is_subband_area_of_interest()
        // we currently use 3.
        const uint filterWidth = 4U;

        var trAr = tilec.resolutions;
        var trPos = 0;
        TcdResolution tr = trAr[trPos], trMax = trAr[numres - 1];

        //Width of the resolution level computed
        var rw = (uint)(tr.x1 - tr.x0);

        //Height of the resolution level computed
        var rh = (uint)(tr.y1 - tr.y0);

        // Compute the intersection of the area of interest, expressed in tile coordinates
        // with the tile coordinates
        var winTcx0 = tilec.win_x0;
        var winTcy0 = tilec.win_y0;
        var winTcx1 = tilec.win_x1;
        var winTcy1 = tilec.win_y1;

        if (trMax.x0 == trMax.x1 || trMax.y0 == trMax.y1)
        {
            return true;
        }

        var sa = SparseArrayInt32.Init(tilec, numres);
        if (sa == null)
            return false;

        if (numres == 1U)
        {
            var ret = sa.read(trMax.win_x0 - (uint)trMax.x0,
                trMax.win_y0 - (uint)trMax.y0,
                trMax.win_x1 - (uint)trMax.x0,
                trMax.win_y1 - (uint)trMax.y0,
                tilec.data_win, 0,
                1, trMax.win_x1 - trMax.win_x0,
                true);
            Debug.Assert(ret);
            return true;
        }
        ulong dataSize = MaxResolution(trAr, (int)numres);
        // overflow check
        // C# 
        if (dataSize > Constants.SizeMax / NbEltsV8 * sizeof(float))
        {
            return false;
        }
        h.Wavelet = new float[dataSize * NbEltsV8];
        v.Wavelet = h.Wavelet;

        for (uint resno = 1; resno < numres; resno++)
        {
            uint i, j;
            /* Window of interest subband-based coordinates */
            uint winLlX0, winLlY0, winLlX1, winLlY1;
            uint winHlX0, winHlX1;
            uint winLhY0, winLhY1;
            /* Window of interest tile-resolution-based coordinates */
            uint winTrX0, winTrX1, winTrY0, winTrY1;
            /* Tile-resolution subband-based coordinates */

            tr = trAr[++trPos];

            h.Sn = (int)rw;
            v.Sn = (int)rh;

            rw = (uint)(tr.x1 - tr.x0);
            rh = (uint)(tr.y1 - tr.y0);

            h.Dn = (int)(rw - (uint)h.Sn);
            h.Cas = tr.x0 % 2;

            v.Dn = (int)(rh - (uint)v.Sn);
            v.Cas = tr.y0 % 2;

            // Get the subband coordinates for the window of interest
            // LL band
            GetBandCoordinates(tilec, resno, 0,
                winTcx0, winTcy0, winTcx1, winTcy1,
                out winLlX0, out winLlY0,
                out winLlX1, out winLlY1);

            // HL band
            GetBandCoordinates(tilec, resno, 1,
                winTcx0, winTcy0, winTcx1, winTcy1,
                out winHlX0, out _, out winHlX1, out _);

            /* LH band */
            GetBandCoordinates(tilec, resno, 2,
                winTcx0, winTcy0, winTcx1, winTcy1,
                out _, out winLhY0, out _, out winLhY1);

            /* Beware: band index for non-LL0 resolution are 0=HL, 1=LH and 2=HH */
            var trLlX0 = (uint)tr.bands[1].x0;
            var trLlY0 = (uint)tr.bands[0].y0;
            var trHlX0 = (uint)tr.bands[0].x0;
            var trLhY0 = (uint)tr.bands[1].y0;

            /* Subtract the origin of the bands for this tile, to the subwindow */
            /* of interest band coordinates, so as to get them relative to the */
            /* tile */
            winLlX0 = MyMath.uint_subs(winLlX0, trLlX0);
            winLlY0 = MyMath.uint_subs(winLlY0, trLlY0);
            winLlX1 = MyMath.uint_subs(winLlX1, trLlX0);
            winLlY1 = MyMath.uint_subs(winLlY1, trLlY0);
            winHlX0 = MyMath.uint_subs(winHlX0, trHlX0);
            winHlX1 = MyMath.uint_subs(winHlX1, trHlX0);
            winLhY0 = MyMath.uint_subs(winLhY0, trLhY0);
            winLhY1 = MyMath.uint_subs(winLhY1, trLhY0);

            SegmentGrow(filterWidth, (uint)h.Sn, ref winLlX0, ref winLlX1);
            SegmentGrow(filterWidth, (uint)h.Dn, ref winHlX0, ref winHlX1);

            SegmentGrow(filterWidth, (uint)v.Sn, ref winLlY0, ref winLlY1);
            SegmentGrow(filterWidth, (uint)v.Dn, ref winLhY0, ref winLhY1);

            /* Compute the tile-resolution-based coordinates for the window of interest */
            if (h.Cas == 0)
            {
                winTrX0 = Math.Min(2 * winLlX0, 2 * winHlX0 + 1);
                winTrX1 = Math.Min(Math.Max(2 * winLlX1, 2 * winHlX1 + 1), rw);
            }
            else
            {
                winTrX0 = Math.Min(2 * winHlX0, 2 * winLlX0 + 1);
                winTrX1 = Math.Min(Math.Max(2 * winHlX1, 2 * winLlX1 + 1), rw);
            }

            if (v.Cas == 0)
            {
                winTrY0 = Math.Min(2 * winLlY0, 2 * winLhY0 + 1);
                winTrY1 = Math.Min(Math.Max(2 * winLlY1, 2 * winLhY1 + 1), rh);
            }
            else
            {
                winTrY0 = Math.Min(2 * winLhY0, 2 * winLlY0 + 1);
                winTrY1 = Math.Min(Math.Max(2 * winLhY1, 2 * winLlY1 + 1), rh);
            }

            h.WinLX0 = winLlX0;
            h.WinLX1 = winLlX1;
            h.WinHX0 = winHlX0;
            h.WinHX1 = winHlX1;
            for (j = 0; j + (NbEltsV8 - 1) < rh; j += NbEltsV8)
            {
                if ((j + (NbEltsV8 - 1) >= winLlY0 && j < winLlY1) ||
                    (j + (NbEltsV8 - 1) >= winLhY0 + (uint)v.Sn &&
                     j < winLhY1 + (uint)v.Sn))
                {
                    v8dwt_interleave_partial_h(h, sa, j, Math.Min(NbEltsV8, rh - j));
                    v8dwt_decode(h);
                    if (!sa.write(winTrX0, j,
                            winTrX1, j + NbEltsV8,
                            h.Wavelet, (int)winTrX0 * NbEltsV8, //C# Wavlet indexing, therefore mul with NB_ELTS_V8
                            NbEltsV8, 1, true))
                    {
                        return false;
                    }
                }
            }

            if (j < rh &&
                ((j + (NbEltsV8 - 1) >= winLlY0 && j < winLlY1) ||
                 (j + (NbEltsV8 - 1) >= winLhY0 + (uint)v.Sn &&
                  j < winLhY1 + (uint)v.Sn)))
            {
                v8dwt_interleave_partial_h(h, sa, j, rh - j);
                v8dwt_decode(h);
                if (!sa.write(winTrX0, j,
                        winTrX1, rh,
                        h.Wavelet, (int)winTrX0 * NbEltsV8,
                        NbEltsV8, 1, true))
                {
                    return false;
                }
            }

            v.WinLX0 = winLlY0;
            v.WinLX1 = winLlY1;
            v.WinHX0 = winLhY0;
            v.WinHX1 = winLhY1;
            for (j = winTrX0; j < winTrX1; j += NbEltsV8)
            {
                var nbElts = Math.Min(NbEltsV8, winTrX1 - j);

                v8dwt_interleave_partial_v(v, sa, j, nbElts);
                v8dwt_decode(v);
                if (!sa.write(j, winTrY0,
                        j + nbElts, winTrY1,
                        v.Wavelet, (int)winTrY0 * NbEltsV8,
                        1, NbEltsV8, true))
                {
                    return false;
                }
            }
        }

        {
            var ret = sa.read(
                trMax.win_x0 - (uint)trMax.x0,
                trMax.win_y0 - (uint)trMax.y0,
                trMax.win_x1 - (uint)trMax.x0,
                trMax.win_y1 - (uint)trMax.y0,
                tilec.data_win, 0,
                1, trMax.win_x1 - trMax.win_x0,
                true);
            Debug.Assert(ret);
        }

        return true;
    }

    /// <summary>
    /// Inverse 9-7 wavelet transform in 2-D.
    /// </summary>
    /// <remarks>2.5.1 - opj_dwt_decode_tile</remarks>
    internal static bool decode_tile_97(TcdTilecomp tilec, uint numres, bool disableMultiThreading)
    {
        var h = new V4dwtLocal();
        var v = new V4dwtLocal();

        var resAr = tilec.resolutions;
        var resArPos = 0;
        var res = resAr[resArPos];

        //Width of the resolution level computed
        var rw = (uint)(res.x1 - res.x0);

        //Height of the resolution level computed
        var rh = (uint)(res.y1 - res.y0);

        var w = tilec.resolutions[tilec.minimum_num_resolutions - 1].x1 
                - tilec.resolutions[tilec.minimum_num_resolutions - 1].x0;

        int numThreads;

        // Not entirely sure for the return code of w == 0 which is triggered per
        // https://github.com/uclouvain/openjpeg/issues/1505
        if (numres == 1U || w == 0)
        {
            return true;
        }
        ThreadPool.GetAvailableThreads(out numThreads, out _);
        numThreads = Math.Min(Environment.ProcessorCount, numThreads);
        if (disableMultiThreading)
            numThreads = 0;

        ulong dataSize = MaxResolution(resAr, (int)numres);
        /* overflow check */
        if (dataSize > Constants.SizeMax / NbEltsV8 * sizeof(float))
        {
            return false;
        }

        h.Wavelet = new float[dataSize * NbEltsV8];
        v.Wavelet = h.Wavelet;

        using (var reset = new ManualResetEvent(false))
        {
            while (--numres != 0)
            {
                var ajAr = tilec.data;
                int aj = 0, j;

                h.Sn = (int)rw;
                v.Sn = (int)rh;

                res = resAr[++resArPos];

                rw = (uint)(res.x1 - res.x0);
                rh = (uint)(res.y1 - res.y0);

                h.Dn = (int)(rw - (uint)h.Sn);
                h.Cas = res.x0 % 2;

                h.WinLX0 = 0;
                h.WinLX1 = (uint)h.Sn;
                h.WinHX0 = 0;
                h.WinHX1 = (uint)h.Dn;

                //C# impl. note
                //Used to convert from float to the "int value" of the
                //float. This is since the org impl. use a int array
                //to store both float and int values
                var fi = new IntOrFloat();

                if (numThreads <= 1 || rh < 2 * NbEltsV8)
                {
                    for (j = 0; j + (NbEltsV8 - 1) < rh; j += NbEltsV8)
                    {
                        v8dwt_interleave_h(h, ajAr, aj, w, NbEltsV8);
                        v8dwt_decode(h);

                        //Copies that back into the aj array.
                        // C# I'm unsure why it's split into two loops, that is probably an
                        // optimalization.
                        for (var k = 0; k < rw; k++)
                        {
                            //C# note: Org. impl stores the wavlet as a struct with eight
                            //floating points. Here it's stored as a continuous array, this
                            //is why we have to multiply k with 8. 
                            var kWavelet = k * NbEltsV8;
                            fi.F = h.Wavelet[kWavelet + 0];
                            ajAr[aj + k] = fi.I;
                            fi.F = h.Wavelet[kWavelet + 1];
                            ajAr[aj + k + w] = fi.I;
                            fi.F = h.Wavelet[kWavelet + 2];
                            ajAr[aj + k + w * 2] = fi.I;
                            fi.F = h.Wavelet[kWavelet + 3];
                            ajAr[aj + k + w * 3] = fi.I;
                        }
                        for (var k = 0; k < rw; k++)
                        {
                            var kWavelet = k * NbEltsV8;
                            fi.F = h.Wavelet[kWavelet + 4];
                            ajAr[aj + k + w * 4] = fi.I;
                            fi.F = h.Wavelet[kWavelet + 5];
                            ajAr[aj + k + w * 5] = fi.I;
                            fi.F = h.Wavelet[kWavelet + 6];
                            ajAr[aj + k + w * 6] = fi.I;
                            fi.F = h.Wavelet[kWavelet + 7];
                            ajAr[aj + k + w * 7] = fi.I;
                        }

                        aj += w * NbEltsV8;
                    }
                }
                else
                {
                    var numJobs = numThreads;

                    if (rh / NbEltsV8 < numJobs)
                    {
                        numJobs = (int)rh / NbEltsV8;
                    }

                    var stepJ = rh / (uint)numJobs / NbEltsV8 * NbEltsV8;

                    reset.Reset();
                    //Alternativly, we can set this to num_jobs and remove the Interlocked.Increment
                    //and the Interlocked.Decrement after the for loop
                    var nThreadWorkers = 1;

                    for (j = 0; j < numJobs; j++)
                    {
                        var job = new Dwt97DecodeHJob(
                            h.Clone(), rw, (uint)w, ajAr, aj,
                            j + 1 == numJobs ? (rh & unchecked((uint)~(NbEltsV8 - 1))) - (uint)j * stepJ : stepJ
                        )
                        {
                            H =
                            {
                                Wavelet = new float[h.Wavelet.Length]
                            }
                        };

                        aj += (int)(w * job.NbRows);

                        Interlocked.Increment(ref nThreadWorkers);
                        ThreadPool.QueueUserWorkItem((x) =>
                        {
                            try { dwt97_decode_h_func((Dwt97DecodeHJob)x); }
                            finally
                            {
                                if (Interlocked.Decrement(ref nThreadWorkers) == 0)
                                    reset.Set();
                            }
                        }, job);
                    }
                    if (Interlocked.Decrement(ref nThreadWorkers) == 0)
                        reset.Set();
                    reset.WaitOne();
                    j = (int)(rh & unchecked((uint)~(NbEltsV8 - 1)));
                }

                if (j < rh)
                {
                        
                    v8dwt_interleave_h(h, ajAr, aj, w, rh - (uint)j);
                    v8dwt_decode(h);

                    for (var k = 0; k < rw; k++)
                    {
                        int kWavelet = k * NbEltsV8, ajk = aj + k;
                        for (uint l = 0; l < rh - j; l++)
                        {
                            fi.F = h.Wavelet[kWavelet + l];
                            ajAr[ajk + w * l] = fi.I;
                        }
                    }
                }

                v.Dn = (int)(rh - (uint)v.Sn);
                v.Cas = res.y0 % 2;
                v.WinLX0 = 0;
                v.WinLX1 = (uint)v.Sn;
                v.WinHX0 = 0;
                v.WinHX1 = (uint)v.Dn;

                aj = 0;
                if (numThreads <= 1 || rw < 2 * NbEltsV8)
                {
                    for (j = (int)rw; j > NbEltsV8 - 1; j -= NbEltsV8)
                    {
                        IntOrFloat faa;
                        faa.I = 0;

                        v8dwt_interleave_v(v, ajAr, aj, w, NbEltsV8);
                        v8dwt_decode(v);
                        for (var k = 0; k < rh; ++k)
                        {
                            Buffer.BlockCopy(v.Wavelet, k * NbEltsV8 * sizeof(float), ajAr, (aj + k * w) * sizeof(float), NbEltsV8 * sizeof(float));
                        }

                        aj += NbEltsV8;
                    }
                }
                else
                {
                    /* "bench_dwt -I" shows that scaling is poor, likely due to RAM
                        transfer being the limiting factor. So limit the number of
                        threads.
                        C# note: I've not run this benchmark
                     */
                    var numJobs = Math.Max(numThreads / 2, 2);
                    //num_jobs = 1;

                    if (rw / NbEltsV8 < numJobs)
                    {
                        numJobs = (int)rw / NbEltsV8;
                    }

                    var stepJ = rw / (uint)numJobs / NbEltsV8 * NbEltsV8;

                    reset.Reset();
                    //Alternativly, we can set this to num_jobs and remove the Interlocked.Increment
                    //and the Interlocked.Decrement after the for loop
                    var nThreadWorkers = 1;

                    for (j = 0; j < numJobs; j++)
                    {
                        var job = new Dwt97DecodeVJob(
                            v.Clone(), rh, (uint)w, ajAr, aj,
                            j + 1 == numJobs ? (rw & unchecked((uint)~(NbEltsV8 - 1))) - (uint)j * stepJ : stepJ
                        )
                        {
                            V =
                            {
                                Wavelet = new float[v.Wavelet.Length]
                            }
                        };

                        aj += (int)job.NbColumns;

                        Interlocked.Increment(ref nThreadWorkers);
                        ThreadPool.QueueUserWorkItem((x) =>
                        {
                            try { dwt97_decode_v_func((Dwt97DecodeVJob)x); }
                            finally
                            {
                                if (Interlocked.Decrement(ref nThreadWorkers) == 0)
                                    reset.Set();
                            }
                        }, job);
                    }
                    if (Interlocked.Decrement(ref nThreadWorkers) == 0)
                        reset.Set();
                    reset.WaitOne();
                }

                //Makes sure not not overflow the array by copying "less than 4 floats".
                if ((rw & (NbEltsV8 - 1)) != 0)
                {
                    j = (int)(rw & (NbEltsV8 - 1));
                    v8dwt_interleave_v(v, ajAr, aj, w, j);
                    v8dwt_decode(v);

                    for (var k = 0; k < rh; ++k)
                        Buffer.BlockCopy(v.Wavelet, k * NbEltsV8 * sizeof(float), ajAr, (aj + k * w) * sizeof(float), j * sizeof(float));
                }
            }
        }

        //{
        //    //C# Debug code to dump out the wavelets. Reason for doing this was to get the same floating point
        //    //   precision as the original impl. 
        //    IntOrFloat faa;
        //    faa.I = 0;

        //    using (var file = new System.IO.StreamWriter("c:/temp/cs_wavelet.txt", append: true))
        //    {
        //        file.Write("Raw float values for {0} wavelets\n", h.win_h_x1 - h.win_h_x0);
        //        for (uint wave_nr = h.win_h_x0; wave_nr < h.win_h_x1; wave_nr++)
        //        {
        //            file.Write("  --- wave {0}\n", wave_nr);
        //            for (int wave_comp = 0; wave_comp < 8; wave_comp++)
        //            {
        //                faa.F = h.wavelet[wave_nr * 8 + wave_comp];
        //                //FP is rounded differently, {1:0.000000}
        //                file.Write(string.Format(System.Globalization.CultureInfo.InvariantCulture, "{1}: {0}\n", faa.I, wave_comp + 1));
        //            }
        //        }
        //    }
        //}

        return true;
    }

    /// <summary>
    /// Inverse 9-7 wavelet transform in 1-D.
    /// </summary>
    /// <remarks>
    /// 2.5
    /// </remarks>
    private static void v8dwt_decode(V4dwtLocal dwt)
    {
        /* BUG_WEIRD_TWO_INVK (look for this identifier in tcd.c) */
        /* Historic value for 2 / opj_invK */
        /* Normally, we should use invK, but if we do so, we have failures in the */
        /* conformance test, due to MSE and peak errors significantly higher than */
        /* accepted value */
        /* Due to using two_invK instead of invK, we have to compensate in tcd.c */
        /* the computation of the stepsize for the non LL subbands */
        const float twoInvK = 1.625732422f;

        int a, b;
        if (dwt.Cas == 0)
        {
            if (!(dwt.Dn > 0 || dwt.Sn > 1))
                return;
            a = 0;
            b = 1;
        }
        else
        {
            if (!(dwt.Sn > 0 || dwt.Dn > 1))
                return;
            a = 1;
            b = 0;
        }

        //C# Snip SSE code

        //C# Both a and b index dwt.wavelet, hench why we have to mul with NB_ELTS_V8
        v8dwt_decode_step1(dwt.Wavelet, a * NbEltsV8, (int)dwt.WinLX0, (int)dwt.WinLX1, K);
        v8dwt_decode_step1(dwt.Wavelet, b * NbEltsV8, (int)dwt.WinHX0, (int)dwt.WinHX1, twoInvK);
        v8dwt_decode_step2(dwt.Wavelet, b * NbEltsV8, (a + 1) * NbEltsV8, (int)dwt.WinLX0, (int)dwt.WinLX1,
            Math.Min(dwt.Sn, dwt.Dn - a), -DwtDelta);
        v8dwt_decode_step2(dwt.Wavelet, a * NbEltsV8, (b + 1) * NbEltsV8, (int)dwt.WinHX0, (int)dwt.WinHX1,
            Math.Min(dwt.Dn, dwt.Sn - b), -DwtGamma);
        v8dwt_decode_step2(dwt.Wavelet, b * NbEltsV8, (a + 1) * NbEltsV8, (int)dwt.WinLX0, (int)dwt.WinLX1,
            Math.Min(dwt.Sn, dwt.Dn - a), -DwtBeta);
        v8dwt_decode_step2(dwt.Wavelet, a * NbEltsV8, (b + 1) * NbEltsV8, (int)dwt.WinHX0, (int)dwt.WinHX1,
            Math.Min(dwt.Dn, dwt.Sn - b), -DwtAlpha);
    }

    /// <summary>
    /// Wavelet decode step 1
    /// </summary>
    /// <remarks>
    /// 2.5 - opj_v8dwt_decode_step1
    /// </remarks>
    /// <param name="wavelet">The array with wavelet</param>
    /// <param name="start">Position of the wavelet in the array</param>
    /// <param name="end">Wavelet count</param>
    /// <param name="c">Some constant</param>
    private static void v8dwt_decode_step1(float[] f, int fw, int start, int end, float c)
    {
        // To be adapted if NB_ELTS_V8 changes
        for (var i = start; i < end; ++i)
        {
            f[fw + i * 2 * 8 + 0] = f[fw + i * 2 * 8 + 0] * c;
            f[fw + i * 2 * 8 + 1] = f[fw + i * 2 * 8 + 1] * c;
            f[fw + i * 2 * 8 + 2] = f[fw + i * 2 * 8 + 2] * c;
            f[fw + i * 2 * 8 + 3] = f[fw + i * 2 * 8 + 3] * c;
            f[fw + i * 2 * 8 + 4] = f[fw + i * 2 * 8 + 4] * c;
            f[fw + i * 2 * 8 + 5] = f[fw + i * 2 * 8 + 5] * c;
            f[fw + i * 2 * 8 + 6] = f[fw + i * 2 * 8 + 6] * c;
            f[fw + i * 2 * 8 + 7] = f[fw + i * 2 * 8 + 7] * c;
        }
    }

    //2.5 - opj_v8dwt_encode_step1
    private static void v8dwt_encode_step1(int[] f, int fw, uint end, float cst)
    {
        var fwi = new IntOrFloat();
        for (var i = 0; i < end; ++i)
        {
            for (var c = 0; c < NbEltsV8; c++)
            {
                fwi.I = f[fw + i * 2 * NbEltsV8 + c];
                fwi.F *= cst;
                f[fw + i * 2 * NbEltsV8 + c] = fwi.I;
            }
        }
    }

    //2.5 - opj_v8dwt_encode_step2
    private static void v8dwt_encode_step2(int[] f, int fl, int fw, uint end, uint m, float cst)
    {
        var imax = Math.Min(end, m);
        var fli = new IntOrFloat();
        var fwi = new IntOrFloat();
        //Snip SSE code

        if (imax > 0)
        {
            for(var c = 0; c < NbEltsV8; c++)
            {
                fli.I = f[fl +  0 * NbEltsV8 + c];
                fwi.I = f[fw +  0 * NbEltsV8 + c];
#if TEST_MATH_MODE
                fli.F = (fli.F + fwi.F) * cst;
#else
                    fli.F = (fli.F + fwi.F) * cst;
#endif
                fwi.I = f[fw + -1 * NbEltsV8 + c];
                fwi.F += fli.F;
                f[fw + -1 * NbEltsV8 + c] = fwi.I;
            }
            fw += 2 * NbEltsV8;
            for(var i = 1; i < imax; i++)
            {
                for (var c = 0; c < NbEltsV8; c++)
                {
                    fli.I = f[fw + -2 * NbEltsV8 + c];
                    fwi.I = f[fw +  0 * NbEltsV8 + c];
#if TEST_MATH_MODE
                    fli.F = (fli.F + fwi.F) * cst;
#else
                        fli.F = (fli.F + fwi.F) * cst;
#endif
                    fwi.I = f[fw + -1 * NbEltsV8 + c];
                    fwi.F += fli.F;
                    f[fw + -1 * NbEltsV8 + c] = fwi.I;
                }
                fw += 2 * NbEltsV8;
            }
        }
        if (m < end)
        {
            Debug.Assert(m + 1 == end);
            for (var c = 0; c < NbEltsV8; c++)
            {
                fwi.I = f[fw + -2 * NbEltsV8 + c];
#if TEST_MATH_MODE
                fli.F = 2 * fwi.F * cst;
#else
                    fli.F = (2 * fwi.F) * cst;
#endif
                fwi.I = f[fw + -1 * NbEltsV8 + c];
                fwi.F += fli.F;
                f[fw + -1 * NbEltsV8 + c] = fwi.I;
            }
        }
    }

    //2.5 - opj_v8dwt_decode_step2
    private static void v8dwt_decode_step2(float[] f, int fl, int fw, int start, int end, int m, float c)
    {
        var imax = Math.Min(end, m);
        if (start > 0)
        {
            fw += 2 * NbEltsV8 * start;
            fl = fw - 2 * NbEltsV8;
        }
        /* To be adapted if NB_ELTS_V8 changes */
        for (var i = start; i < imax; ++i)
        {
            //if (i == 5)
            //{
            //    IntOrFloat faa;
            //    faa.I = 0;

            //    faa.F = f[fw - 7];
            //    faa.F = f[fl + 1];
            //    faa.F = f[fw + 1];
            //    faa.F = c;

            //    float prob_f = (float)(f[fw - 7] + (float)((float)(f[fl + 1] + f[fw + 1]) * (float)c));
            //    faa.F = prob_f;

            //    i = i;
            //}

#if TEST_MATH_MODE
            //C# We have a problem. AFAICT, C# likes to do double precision math
            //   Liberal use of (float) seems to fix the issue. 
            f[fw - 8] = f[fw - 8] + (f[fl + 0] + f[fw + 0]) * c;
            f[fw - 7] = f[fw - 7] + (f[fl + 1] + f[fw + 1]) * c;
            f[fw - 6] = f[fw - 6] + (f[fl + 2] + f[fw + 2]) * c;
            f[fw - 5] = f[fw - 5] + (f[fl + 3] + f[fw + 3]) * c;
            f[fw - 4] = f[fw - 4] + (f[fl + 4] + f[fw + 4]) * c;
            f[fw - 3] = f[fw - 3] + (f[fl + 5] + f[fw + 5]) * c;
            f[fw - 2] = f[fw - 2] + (f[fl + 6] + f[fw + 6]) * c;
            f[fw - 1] = f[fw - 1] + (f[fl + 7] + f[fw + 7]) * c; 
#else
                //It's not a bug to have greater precision, it's a problem for verifying
                //the result against original libary.
                f[fw - 8] = f[fw - 8] + ((f[fl + 0] + f[fw + 0]) * c);
                f[fw - 7] = f[fw - 7] + ((f[fl + 1] + f[fw + 1]) * c);
                f[fw - 6] = f[fw - 6] + ((f[fl + 2] + f[fw + 2]) * c);
                f[fw - 5] = f[fw - 5] + ((f[fl + 3] + f[fw + 3]) * c);
                f[fw - 4] = f[fw - 4] + ((f[fl + 4] + f[fw + 4]) * c);
                f[fw - 3] = f[fw - 3] + ((f[fl + 5] + f[fw + 5]) * c);
                f[fw - 2] = f[fw - 2] + ((f[fl + 6] + f[fw + 6]) * c);
                f[fw - 1] = f[fw - 1] + ((f[fl + 7] + f[fw + 7]) * c);
#endif
            fl = fw;
            fw += 2 * NbEltsV8;
        }
        if (m < end)
        {
            Debug.Assert(m + 1 == end);
            c += c;
#if TEST_MATH_MODE
            f[fw - 8] = f[fw - 8] + f[fl + 0] * c;
            f[fw - 7] = f[fw - 7] + f[fl + 1] * c;
            f[fw - 6] = f[fw - 6] + f[fl + 2] * c;
            f[fw - 5] = f[fw - 5] + f[fl + 3] * c;
            f[fw - 4] = f[fw - 4] + f[fl + 4] * c;
            f[fw - 3] = f[fw - 3] + f[fl + 5] * c;
            f[fw - 2] = f[fw - 2] + f[fl + 6] * c;
            f[fw - 1] = f[fw - 1] + f[fl + 7] * c;
#else
                f[fw - 8] = f[fw - 8] + (f[fl + 0] * c);
                f[fw - 7] = f[fw - 7] + (f[fl + 1] * c);
                f[fw - 6] = f[fw - 6] + (f[fl + 2] * c);
                f[fw - 5] = f[fw - 5] + (f[fl + 3] * c);
                f[fw - 4] = f[fw - 4] + (f[fl + 4] * c);
                f[fw - 3] = f[fw - 3] + (f[fl + 5] * c);
                f[fw - 2] = f[fw - 2] + (f[fl + 6] * c);
                f[fw - 1] = f[fw - 1] + (f[fl + 7] * c);
#endif
        }
    }

    //2.5
    private static void EncodeStep1Combined(int[] f, int fw, uint itersC1, uint itersC2, float c1, float c2)
    {
        var fw1 = new IntOrFloat();
        uint i = 0;
        var itersCommon = Math.Min(itersC1, itersC2);
        Debug.Assert(Math.Abs((int)itersC1 - (int)itersC2) <= 1);
        for(; i + 3 < itersCommon; i += 4)
        {
            fw1.I = f[fw + 0];
            fw1.F *= c1;
            f[fw + 0] = fw1.I;

            fw1.I = f[fw + 1];
            fw1.F *= c2;
            f[fw + 1] = fw1.I;

            fw1.I = f[fw + 2];
            fw1.F *= c1;
            f[fw + 2] = fw1.I;

            fw1.I = f[fw + 3];
            fw1.F *= c2;
            f[fw + 3] = fw1.I;

            fw1.I = f[fw + 4];
            fw1.F *= c1;
            f[fw + 4] = fw1.I;

            fw1.I = f[fw + 5];
            fw1.F *= c2;
            f[fw + 5] = fw1.I;

            fw1.I = f[fw + 6];
            fw1.F *= c1;
            f[fw + 6] = fw1.I;

            fw1.I = f[fw + 7];
            fw1.F *= c2;
            f[fw + 7] = fw1.I;

            fw += 8;
        }
        for(; i < itersCommon; i++)
        {
            fw1.I = f[fw + 0];
            fw1.F *= c1;
            f[fw + 0] = fw1.I;

            fw1.I = f[fw + 1];
            fw1.F *= c2;
            f[fw + 1] = fw1.I;

            fw += 2;
        }
        if (i < itersC1)
        {
            fw1.I = f[fw + 0];
            fw1.F *= c1;
            f[fw + 0] = fw1.I;
        }
        else if (i < itersC2)
        {
            fw1.I = f[fw + 1];
            fw1.F *= c2;
            f[fw + 1] = fw1.I;
        }
    }

    //2.5
    private static void EncodeStep2(int[] f, int fl, int fw, uint end, uint m, float c)
    {
        IntOrFloat fw1 = new IntOrFloat(), fw2 = new IntOrFloat();
        var imax = Math.Min(end, m);
        if (imax > 0)
        {
            //fw[-1] += (fl[0] + fw[0]) * c;
            fw1.I = f[fl + 0];
            fw2.I = f[fw + 0];
#if TEST_MATH_MODE
            fw2.F = (fw1.F + fw2.F) * c;
#else
                fw2.F = (fw1.F + fw2.F) * c;
#endif
            fw1.I = f[fw - 1];
            fw1.F += fw2.F;
            f[fw - 1] = fw1.I;
            fw += 2;
            var i = 1;
            for(; i + 3 < imax; i += 4)
            {
                //fw[-1] += (fw[-2] + fw[0]) * c;
                fw1.I = f[fw - 2];
                fw2.I = f[fw + 0];
#if TEST_MATH_MODE
                fw2.F = (fw1.F + fw2.F) * c;
#else
                    fw2.F = (fw1.F + fw2.F) * c;
#endif
                fw1.I = f[fw - 1];
                fw1.F += fw2.F;
                f[fw - 1] = fw1.I;

                //fw[1] += (fw[0] + fw[2]) * c;
                fw1.I = f[fw + 0];
                fw2.I = f[fw + 2];
#if TEST_MATH_MODE
                fw2.F = (fw1.F + fw2.F) * c;
#else
                    fw2.F = (fw1.F + fw2.F) * c;
#endif
                fw1.I = f[fw + 1];
                fw1.F += fw2.F;
                f[fw + 1] = fw1.I;

                //fw[3] += (fw[2] + fw[4]) * c;
                fw1.I = f[fw + 2];
                fw2.I = f[fw + 4];
#if TEST_MATH_MODE
                fw2.F = (fw1.F + fw2.F) * c;
#else
                    fw2.F = (fw1.F + fw2.F) * c;
#endif
                fw1.I = f[fw + 3];
                fw1.F += fw2.F;
                f[fw + 3] = fw1.I;

                //fw[5] += (fw[4] + fw[6]) * c;
                fw1.I = f[fw + 4];
                fw2.I = f[fw + 6];
#if TEST_MATH_MODE
                fw2.F = (fw1.F + fw2.F) * c;
#else
                    fw2.F = (fw1.F + fw2.F) * c;
#endif
                fw1.I = f[fw + 5];
                fw1.F += fw2.F;
                f[fw + 5] = fw1.I;

                fw += 8;
            }
            for(; i < imax; i++)
            {
                //fw[-1] += (fw[-2] + fw[0]) * c
                fw1.I = f[fw - 2];
                fw2.I = f[fw + 0];
#if TEST_MATH_MODE
                fw2.F = (fw1.F + fw2.F) * c;
#else
                    fw2.F = (fw1.F + fw2.F) * c;
#endif
                fw1.I = f[fw - 1];
                fw1.F += fw2.F;
                f[fw - 1] = fw1.I;

                fw += 2;
            }
        }
        if (m < end)
        {
            Debug.Assert(m + 1 == end);
            //fw[-1] += (2 * fw[-2]) * c;
            fw2.I = f[fw - 2];
#if TEST_MATH_MODE
            fw2.F = 2 * fw2.F * c;
#else
                fw2.F = (2 * fw2.F) * c;
#endif
            fw1.I = f[fw - 1];
            fw1.F += fw2.F;
            f[fw - 1] = fw1.I;
        }
    }

    //2.5
    private static bool decode_tile(TcdTilecomp tilec, uint numres, bool disableMultiThreading)
    {
        var h = new DwtLocal();
        var v = new DwtLocal();

        var trAr = tilec.resolutions;
        var trPos = 0;
        TcdResolution tr = trAr[trPos], trMax = trAr[numres - 1];

        //Width of the resolution level computed
        var rw = (uint)(tr.x1 - tr.x0);

        //Height of the resolution level computed
        var rh = (uint)(tr.y1 - tr.y0);

        var w = (uint)(tilec.resolutions[tilec.minimum_num_resolutions -
                                         1].x1 -
                       tilec.resolutions[tilec.minimum_num_resolutions - 1].x0);

        int numThreads;

        if (numres == 1U)
        {
            return true;
        }
        ThreadPool.GetAvailableThreads(out numThreads, out _);
        numThreads = disableMultiThreading ? 1 : Math.Min(Environment.ProcessorCount, numThreads);

        ulong hMemSize = MaxResolution(trAr, (int)numres);
        /* overflow check */
        if (hMemSize > Constants.SizeMax / ParallelCols53 / sizeof(int))
        {
            return false;
        }
        // We need PARALLEL_COLS_53 times the height of the array,
        // since for the vertical pass
        // we process PARALLEL_COLS_53 columns at a time
        hMemSize *= ParallelCols53;
        h.Mem = new int[hMemSize];
        v.Mem = h.Mem;

        using (var reset = new ManualResetEvent(false))
        {
            while (--numres != 0)
            {
                var tiledp = 0;
                uint j;

                tr = trAr[++trPos];
                h.Sn = (int)rw;
                v.Sn = (int)rh;

                rw = (uint)(tr.x1 - tr.x0);
                rh = (uint)(tr.y1 - tr.y0);

                h.Dn = (int)(rw - (uint)h.Sn);
                h.Cas = tr.x0 % 2;

                if (numThreads <= 1 || rh <= 1)
                {
                    for (j = 0; j < rh; ++j)
                    {
                        idwt53_h(h, tilec.data, tiledp + (int)(j * w));
                    }
                }
                else
                {
                    var numJobs = numThreads;

                    if (rh < numJobs)
                    {
                        numJobs = (int)rh;
                    }

                    var stepJ = rh / (uint)numJobs;
                        
                    reset.Reset();
                    //Alternativly, we can set this to num_jobs and remove the Interlocked.Increment
                    //and the Interlocked.Decrement after the for loop
                    var nThreadWorkers = 1;

                    for (j = 0; j < numJobs; j++)
                    {
                        var maxJ = (j + 1U) * stepJ; // this will overflow
                        if (j == numJobs - 1) //So we clamp max_j
                            maxJ = rh;
                        var job = new DecodeHJob(
                            h.Clone(), rw, w, tilec.data, tiledp,
                            j * stepJ, maxJ
                        )
                        {
                            H =
                            {
                                Mem = new int[hMemSize]
                            }
                        };

                        Interlocked.Increment(ref nThreadWorkers);
                        ThreadPool.QueueUserWorkItem((x) =>
                        {
                            try { decode_h_func((DecodeHJob)x); }
                            finally
                            {
                                if (Interlocked.Decrement(ref nThreadWorkers) == 0)
                                    reset.Set();
                            }
                        }, job);
                    }
                    if (Interlocked.Decrement(ref nThreadWorkers) == 0)
                        reset.Set();
                    reset.WaitOne();
                }

                v.Dn = (int)(rh - (uint)v.Sn);
                v.Cas = tr.y0 % 2;
 
                if (numThreads <= 1 || rw <= 1)
                {
                    for (j = 0; j + ParallelCols53 <= rw; j += ParallelCols53)
                    {
                        idwt53_v(v, tilec.data, tiledp + (int)j, (int)w, (int)ParallelCols53);
                    }
                    if (j < rw)
                        idwt53_v(v, tilec.data, tiledp + (int)j, (int)w, (int)(rw - j));
                }
                else
                {
                    var numJobs = numThreads;

                    if (rw < numJobs)
                    {
                        numJobs = (int)rw;
                    }

                    var stepJ = rw / (uint)numJobs;

                    reset.Reset();
                    var nThreadWorkers = 1;

                    for (j = 0; j < numJobs; j++)
                    {
                        var maxJ = (j + 1U) * stepJ; // this can overflow
                        if (j == numJobs - 1)
                            maxJ = rw;
                        var job = new DecodeVJob(
                            v.Clone(), rh, w, tilec.data, tiledp,
                            j * stepJ, maxJ
                        )
                        {
                            V =
                            {
                                Mem = new int[hMemSize]
                            }
                        };

                        Interlocked.Increment(ref nThreadWorkers);
                        ThreadPool.QueueUserWorkItem((x) =>
                        {
                            try { decode_v_func((DecodeVJob)x); }
                            finally
                            {
                                if (Interlocked.Decrement(ref nThreadWorkers) == 0)
                                    reset.Set();
                            }
                        }, job);
                    }
                    if (Interlocked.Decrement(ref nThreadWorkers) == 0)
                        reset.Set();
                    reset.WaitOne();
                }
            }
        }
        return true;
    }

    //2.5
    private static bool DecodePartialTile(TcdTilecomp tilec, uint numres)
    {
        var h = new DwtLocal();
        var v = new DwtLocal();
        // This value matches the maximum left/right extension given in tables
        // F.2 and F.3 of the standard.
        const uint filterWidth = 2U;

        var trAr = tilec.resolutions;
        var trPos = 0;
        TcdResolution tr = trAr[trPos], trMax = trAr[numres - 1];

        //Width of the resolution level computed
        var rw = (uint)(tr.x1 - tr.x0);

        //Height of the resolution level computed
        var rh = (uint)(tr.y1 - tr.y0);

        // Compute the intersection of the area of interest, expressed in tile coordinates
        // with the tile coordinates
        var winTcx0 = tilec.win_x0;
        var winTcy0 = tilec.win_y0;
        var winTcx1 = tilec.win_x1;
        var winTcy1 = tilec.win_y1;

        if (trMax.x0 == trMax.x1 || trMax.y0 == trMax.y1)
        {
            return true;
        }

        var sa = SparseArrayInt32.Init(tilec, numres);
        if (sa == null)
            return false;

        if (numres == 1U)
        {
            var ret = sa.read(trMax.win_x0 - (uint)trMax.x0,
                trMax.win_y0 - (uint)trMax.y0,
                trMax.win_x1 - (uint)trMax.x0,
                trMax.win_y1 - (uint)trMax.y0,
                tilec.data_win, 0,
                1, trMax.win_x1 - trMax.win_x0,
                true);
            Debug.Assert(ret);
            return true;
        }
        ulong hMemSize = MaxResolution(trAr, (int) numres);
        // overflow check
        // in vertical pass, we process 4 columns at a time
        if (hMemSize > Constants.SizeMax / (4 * sizeof(int)))
        {
            return false;
        }
        hMemSize *= 4;
        h.Mem = new int[hMemSize];
        v.Mem = h.Mem;

        for (uint resno = 1; resno < numres; resno++)
        {
            /* Window of interest subband-based coordinates */
            uint winLlX0, winLlY0, winLlX1, winLlY1;
            uint winHlX0, winHlX1;
            uint winLhY0, winLhY1;
            /* Window of interest tile-resolution-based coordinates */
            uint winTrX0, winTrX1, winTrY0, winTrY1;
            /* Tile-resolution subband-based coordinates */

            tr = trAr[++trPos];

            h.Sn = (int)rw;
            v.Sn = (int)rh;

            rw = (uint)(tr.x1 - tr.x0);
            rh = (uint)(tr.y1 - tr.y0);

            h.Dn = (int)(rw - (uint)h.Sn);
            h.Cas = tr.x0 % 2;

            v.Dn = (int)(rh - (uint)v.Sn);
            v.Cas = tr.y0 % 2;

            // Get the subband coordinates for the window of interest
            // LL band
            GetBandCoordinates(tilec, resno, 0,
                winTcx0, winTcy0, winTcx1, winTcy1,
                out winLlX0, out winLlY0,
                out winLlX1, out winLlY1);

            // HL band
            GetBandCoordinates(tilec, resno, 1,
                winTcx0, winTcy0, winTcx1, winTcy1,
                out winHlX0, out _, out winHlX1, out _);

            /* LH band */
            GetBandCoordinates(tilec, resno, 2,
                winTcx0, winTcy0, winTcx1, winTcy1,
                out _, out winLhY0, out _, out winLhY1);

            /* Beware: band index for non-LL0 resolution are 0=HL, 1=LH and 2=HH */
            var trLlX0 = (uint)tr.bands[1].x0;
            var trLlY0 = (uint)tr.bands[0].y0;
            var trHlX0 = (uint)tr.bands[0].x0;
            var trLhY0 = (uint)tr.bands[1].y0;

            /* Subtract the origin of the bands for this tile, to the subwindow */
            /* of interest band coordinates, so as to get them relative to the */
            /* tile */
            winLlX0 = MyMath.uint_subs(winLlX0, trLlX0);
            winLlY0 = MyMath.uint_subs(winLlY0, trLlY0);
            winLlX1 = MyMath.uint_subs(winLlX1, trLlX0);
            winLlY1 = MyMath.uint_subs(winLlY1, trLlY0);
            winHlX0 = MyMath.uint_subs(winHlX0, trHlX0);
            winHlX1 = MyMath.uint_subs(winHlX1, trHlX0);
            winLhY0 = MyMath.uint_subs(winLhY0, trLhY0);
            winLhY1 = MyMath.uint_subs(winLhY1, trLhY0);

            SegmentGrow(filterWidth, (uint)h.Sn, ref winLlX0, ref winLlX1);
            SegmentGrow(filterWidth, (uint)h.Dn, ref winHlX0, ref winHlX1);

            SegmentGrow(filterWidth, (uint)v.Sn, ref winLlY0, ref winLlY1);
            SegmentGrow(filterWidth, (uint)v.Dn, ref winLhY0, ref winLhY1);

            /* Compute the tile-resolution-based coordinates for the window of interest */
            if (h.Cas == 0)
            {
                winTrX0 = Math.Min(2 * winLlX0, 2 * winHlX0 + 1);
                winTrX1 = Math.Min(Math.Max(2 * winLlX1, 2 * winHlX1 + 1), rw);
            }
            else
            {
                winTrX0 = Math.Min(2 * winHlX0, 2 * winLlX0 + 1);
                winTrX1 = Math.Min(Math.Max(2 * winHlX1, 2 * winLlX1 + 1), rw);
            }

            if (v.Cas == 0)
            {
                winTrY0 = Math.Min(2 * winLlY0, 2 * winLhY0 + 1);
                winTrY1 = Math.Min(Math.Max(2 * winLlY1, 2 * winLhY1 + 1), rh);
            }
            else
            {
                winTrY0 = Math.Min(2 * winLhY0, 2 * winLlY0 + 1);
                winTrY1 = Math.Min(Math.Max(2 * winLhY1, 2 * winLlY1 + 1), rh);
            }

            for (uint j = 0; j < rh; ++j)
            {
                if ((j >= winLlY0 && j < winLlY1) ||
                    (j >= winLhY0 + (uint)v.Sn && j < winLhY1 + (uint)v.Sn))
                {
                    // Avoids dwt.c:1584:44 (in opj_dwt_decode_partial_1): runtime error:
                    // signed integer overflow: -1094795586 + -1094795586 cannot be represented in type 'int'
                    // on opj_decompress -i  ../../openjpeg/MAPA.jp2 -o out.tif -d 0,0,256,256
                    // This is less extreme than memsetting the whole buffer to 0
                    // although we could potentially do better with better handling of edge conditions
                    if (winTrX1 >= 1 && winTrX1 < rw)
                    {
                        h.Mem[winTrX1 - 1] = 0;
                    }
                    if (winTrX1 < rw)
                    {
                        h.Mem[winTrX1] = 0;
                    }

                    interleave_partial_h(h.Mem,
                        h.Cas,
                        sa,
                        j,
                        (uint)h.Sn,
                        winLlX0,
                        winLlX1,
                        winHlX0,
                        winHlX1);
                    decode_partial_1(h.Mem, h.Dn, h.Sn, h.Cas,
                        (int)winLlX0,
                        (int)winLlX1,
                        (int)winHlX0,
                        (int)winHlX1);
                    if (!sa.write(winTrX0, j,
                            winTrX1, j + 1,
                            h.Mem, (int)winTrX0,
                            1, 0, true))
                    {
                        return false;
                    }
                }
            }

            for (var i = winTrX0; i < winTrX1;)
            {
                var nbCols = Math.Min(4U, winTrX1 - i);
                interleave_partial_v(v.Mem,
                    v.Cas,
                    sa,
                    i,
                    nbCols,
                    (uint)v.Sn,
                    winLlY0,
                    winLlY1,
                    winLhY0,
                    winLhY1);
                decode_partial_1_parallel(v.Mem, nbCols, v.Dn, v.Sn, v.Cas,
                    (int)winLlY0,
                    (int)winLlY1,
                    (int)winLhY0,
                    (int)winLhY1);
                if (!sa.write(i, winTrY0,
                        i + nbCols, winTrY1,
                        v.Mem, 4 * (int)winTrY0,
                        1, 4, true))
                {
                    return false;
                }

                i += nbCols;
            }
        }

        {
            var ret = sa.read(
                trMax.win_x0 - (uint)trMax.x0,
                trMax.win_y0 - (uint)trMax.y0,
                trMax.win_x1 - (uint)trMax.x0,
                trMax.win_y1 - (uint)trMax.y0,
                tilec.data_win, 0,
                1, trMax.win_x1 - trMax.win_x0,
                true);
            Debug.Assert(ret);
        }

        return true;
    }

    //2.5 - opj_dwt_interleave_partial_h
    private static void interleave_partial_h(
        int[] dest,
        int cas,
        SparseArrayInt32 sa,
        uint saLine,
        uint sn,
        uint winLX0,
        uint winLX1,
        uint winHX0,
        uint winHX1)
    {
        var ret = sa.read(winLX0, saLine,
            winLX1, saLine + 1,
            dest, cas + 2 * (int)winLX0,
            2, 0, true);
        Debug.Assert(ret);
        ret = sa.read(sn + winHX0, saLine,
            sn + winHX1, saLine + 1,
            dest, 1 - cas + 2 * (int)winHX0,
            2, 0, true);
        Debug.Assert(ret);
    }

    //2.5 - opj_dwt_interleave_partial_v
    private static void interleave_partial_v(
        int[] dest,
        int cas,
        SparseArrayInt32 sa,
        uint saCol,
        uint nbCols,
        uint sn,
        uint winLY0,
        uint winLY1,
        uint winHY0,
        uint winHY1)
    {
        var ret = sa.read(saCol, winLY0,
            saCol + nbCols, winLY1,
            dest, cas * 4 + 2 * 4 * (int)winLY0,
            1, 2 * 4, true);
        Debug.Assert(ret);
        ret = sa.read(saCol, sn + winHY0,
            saCol + nbCols, sn + winHY1,
            dest, (1 - cas) * 4 + 2 * 4 * (int)winHY0,
            1, 2 * 4, true);
        Debug.Assert(ret);
    }

    /// <remarks>
    /// 2.5 - opj_dwt_decode_partial_1
    /// 
    /// OPJ_S(i) => a[(i)*2]
    /// OPJ_D(i) => a[(1+(i)*2)]
    /// OPJ_S_(i) => ((i)<0?OPJ_S(0):((i)>=sn?OPJ_S(sn-1):OPJ_S(i)))
    /// OPJ_D_(i) => ((i)<0?OPJ_D(0):((i)>=dn?OPJ_D(dn-1):OPJ_D(i)))
    /// OPJ_SS_(i) => ((i)<0?OPJ_S(0):((i)>=dn?OPJ_S(dn-1):OPJ_S(i)))
    /// OPJ_DD_(i) => ((i)<0?OPJ_D(0):((i)>=sn?OPJ_D(sn-1):OPJ_D(i)))
    /// 
    /// Substituted:
    /// OPJ_D_(i) => ((i) < 0 ? a[1] : ((i) >= dn ? a[(1 + (dn - 1) * 2)] : a[(1 + (i) * 2)]))
    /// OPJ_S_(i) => ((i)<0?a[0]:((i)>=sn?a[(sn-1)*2]:a[(i)*2]))
    /// OPJ_SS_(i) => ((i)<0?a[0]:((i)>=dn?a[(dn-1)*2]:a[(i)*2]))
    /// OPJ_DD_(i) => ((i)<0?a[1]:((i)>=sn?a[(1+(sn-1)*2)]:a[(1+(i)*2)]))
    /// </remarks>
    private static void decode_partial_1(int[] a, int dn, int sn,
        int cas,
        int winLX0,
        int winLX1,
        int winHX0,
        int winHX1)
    {
        int i;

        if (cas == 0)
        {
            if (dn > 0 || sn > 1)
            { 
                i = winLX0;
                if (i < winLX1)
                {
                    /* Left-most case */
                    a[i * 2] -= ((i - 1 < 0 ? a[1 + 0 * 2] : i - 1 >= dn ? a[1 + (dn - 1) * 2] : a[1 + (i - 1) * 2])
                                 + (i < 0 ? a[1 + 0 * 2] : i >= dn ? a[1 + (dn - 1) * 2] : a[1 + i * 2])
                                 + 2) >> 2;
                    i++;

                    var iMax = winLX1;
                    if (iMax > dn)
                    {
                        iMax = dn;
                    }
                    for (; i < iMax; i++)
                    {
                        /* No bound checking */
                        a[i * 2] -= (a[1 + (i - 1) * 2] + a[1 + i * 2] + 2) >> 2;
                    }
                    for (; i < winLX1; i++)
                    {
                        /* Right-most case */
                        a[i * 2] -= ((i - 1 < 0 ? a[1 + 0 * 2] : i - 1 >= dn ? a[1 + (dn - 1) * 2] : a[1 + (i - 1) * 2])
                                     + (i < 0 ? a[1 + 0 * 2] : i >= dn ? a[1 + (dn - 1) * 2] : a[1 + i * 2])
                                     + 2) >> 2;
                    }
                }

                i = winHX0;
                if (i < winHX1)
                {
                    var iMax = winHX1;
                    if (iMax >= sn)
                    {
                        iMax = sn - 1;
                    }
                    for (; i < iMax; i++)
                    {
                        /* No bound checking */
                        a[1 + i * 2] += (a[i * 2] + a[(i + 1) * 2]) >> 1;
                    }
                    for (; i < winHX1; i++)
                    {
                        /* Right-most case */
                        a[1 + i * 2] += ((i < 0 ? a[0] : i >= sn ? a[(sn - 1) * 2] : a[i * 2])
                                         + (i + 1 < 0 ? a[0] : i + 1 >= sn ? a[(sn - 1) * 2] : a[(i + 1) * 2])
                            ) >> 1;
                    }
                }
            }
        }
        else
        {
            if (sn == 0 && dn == 1)
            {        /* NEW :  CASE ONE ELEMENT */
                a[0] /= 2;
            }
            else
            {
                for (i = winLX0; i < winLX1; i++)
                {
                    a[1 + i * 2] = MyMath.int_sub_no_overflow(
                        a[1 + i * 2],
                        MyMath.int_add_no_overflow(
                            MyMath.int_add_no_overflow(
                                i < 0 ? a[0] : i >= dn ? a[(dn - 1) * 2] : a[i * 2],
                                i + 1 < 0 ? a[0] : i + 1 >= dn ? a[(dn - 1) * 2] : a[(i + 1) * 2]),
                            2) >> 2);
                }
                for (i = winHX0; i < winHX1; i++)
                {
                    a[i * 2] = MyMath.int_add_no_overflow(
                        a[i * 2],
                        MyMath.int_add_no_overflow(
                            i < 0 ? a[1] : i >= sn ? a[1 + (sn - 1) * 2] : a[1 + i * 2],
                            i - 1 < 0 ? a[1] : i - 1 >= sn ? a[1 + (sn - 1) * 2] : a[1 + (i - 1) * 2])
                        >> 1);
                }
            }
        }
    }

    /// <remarks>
    /// 2.5
    /// 
    /// Todo: Make this more readable. Look at the OpenJpeg 2.1 C# impl. 
    ///       where the code is much more readable (Decode_1)
    /// 
    /// OPJ_S_off(i,off) => a[(uint)(i)*2*4+off]
    /// OPJ_D_off(i,off) => a[(1+(uint)(i)*2)*4+off]
    /// OPJ_S__off(i,off) => ((i)<0?OPJ_S_off(0,off):((i)>=sn?OPJ_S_off(sn-1,off):OPJ_S_off(i,off)))
    /// OPJ_D__off(i,off) => ((i)<0?OPJ_D_off(0,off):((i)>=dn?OPJ_D_off(dn-1,off):OPJ_D_off(i,off)))
    /// OPJ_SS__off(i,off) => ((i)<0?OPJ_S_off(0,off):((i)>=dn?OPJ_S_off(dn-1,off):OPJ_S_off(i,off)))
    /// OPJ_DD__off(i,off) => ((i)<0?OPJ_D_off(0,off):((i)>=sn?OPJ_D_off(sn-1,off):OPJ_D_off(i,off)))
    /// 
    /// Substituted:
    /// OPJ_S__off(i,off) => ((i)<0?a[off]:((i)>=sn?a[(uint)(sn-1)*2*4+off]:a[(uint)(i)*2*4+off]))
    /// OPJ_D__off(i,off) => ((i)<0?a[1+off]:((i)>=dn?a[(1+(uint)(dn-1)*2)*4+off]:a[(1+(uint)(i)*2)*4+off]))
    /// OPJ_SS__off(i,off) => ((i)<0?a[off]:((i)>=dn?a[(uint)(dn-1)*2*4+off]:a[(uint)(i)*2*4+off]))
    /// OPJ_DD__off(i,off) => ((i)<0?a[1+off]:((i)>=sn?a[(1+(uint)(sn-1)*2)*4+off]:a[(1+(uint)(i)*2)*4+off])
    /// </remarks>
    private static void decode_partial_1_parallel(
        int[] a,
        uint nbCols,
        int dn, int sn,
        int cas,
        int winLX0,
        int winLX1,
        int winHX0,
        int winHX1)
    {
        int i;
        uint off;

        if (cas == 0)
        {
            if (dn > 0 || sn > 1)
            {
                i = winLX0;
                if (i < winLX1)
                {
                    /* Left-most case */
                    for (off = 0; off < 4; off++)
                    {
                        a[(uint)i * 2 * 4 + off] -= (
                            (i - 1 < 0 ? a[1 * 4 + off] : i - 1 >= dn ? a[(1 + (uint)(dn - 1) * 2) * 4 + off] : a[(1 + (uint)(i - 1) * 2) * 4 + off]) 
                            + (i < 0 ? a[1 * 4 + off] : i >= dn ? a[(1 + (uint)(dn - 1) * 2) * 4 + off] : a[(1 + (uint)i * 2) * 4 + off]) 
                            + 2
                        ) >> 2;
                    }
                    i++;

                    var iMax = winLX1;
                    if (iMax > dn)
                    {
                        iMax = dn;
                    }

                    //Snip SSE2 code. 

                    for (; i < iMax; i++)
                    {
                        /* No bound checking */
                        for (off = 0; off < 4; off++)
                        {
                            a[(uint)i * 2 * 4 + off] -= (
                                a[(1 + (uint)(i - 1) * 2) * 4 + off] 
                                + a[(1 + (uint)i * 2) * 4 + off] 
                                + 2
                            ) >> 2;
                        }
                    }
                    for (; i < winLX1; i++)
                    {
                        /* Right-most case */
                        for (off = 0; off < 4; off++)
                        {
                            a[(uint)i * 2 * 4 + off] -= (
                                (i - 1 < 0 ? a[1 * 4 + off] : i - 1 >= dn ? a[(1 + (uint)(dn - 1) * 2) * 4 + off] : a[(1 + (uint)(i - 1) * 2) * 4 + off]) 
                                + (i < 0 ? a[1 * 4 + off] : i >= dn ? a[(1 + (uint)(dn - 1) * 2) * 4 + off] : a[(1 + (uint)i * 2) * 4 + off]) 
                                + 2
                            ) >> 2;
                        }
                    }
                }

                i = winHX0;
                if (i < winHX1)
                {
                    var iMax = winHX1;
                    if (iMax >= sn)
                    {
                        iMax = sn - 1;
                    }

                    //Snip SSE2 code. 

                    for (; i < iMax; i++)
                    {
                        /* No bound checking */
                        for (off = 0; off < 4; off++)
                        {
                            a[(1 + (uint)i * 2) * 4 + off] += (
                                a[(uint)i * 2 * 4 + off]
                                + a[(uint)(i + 1) * 2 * 4 + off]
                            ) >> 1;
                        }
                    }
                    for (; i < winHX1; i++)
                    {
                        /* Right-most case */
                        for (off = 0; off < 4; off++)
                        {
                            a[(1 + (uint)i * 2) * 4 + off] += (
                                (i < 0 ? a[off] : i >= sn ? a[(uint)(sn - 1) * 2 * 4 + off] : a[(uint)i * 2 * 4 + off]) 
                                + (i + 1 < 0 ? a[off] : i + 1 >= sn ? a[(uint)(sn - 1) * 2 * 4 + off] : a[(uint)(i + 1) * 2 * 4 + off])
                            ) >> 1;
                        }
                    }
                }
            }
        }
        else
        {
            if (sn == 0 && dn == 1)
            {        /* NEW :  CASE ONE ELEMENT */
                for (off = 0; off < 4; off++)
                    a[off] /= 2;
            }
            else
            {
                for (i = winLX0; i < winLX1; i++)
                {
                    for (off = 0; off < 4; off++)
                        a[(1 + (uint)i * 2) * 4 + off] = 
                            MyMath.int_sub_no_overflow(
                                a[(1 + (uint)i * 2) * 4 + off],
                                MyMath.int_add_no_overflow(
                                    MyMath.int_add_no_overflow(
                                        i < 0 ? a[(uint)0 * 2 * 4 + off] : i >= dn ? a[(uint)(dn - 1) * 2 * 4 + off] : a[(uint)i * 2 * 4 + off], 
                                        i + 1 < 0 ? a[(uint)0 * 2 * 4 + off] : i + 1 >= dn ? a[(uint)(dn - 1) * 2 * 4 + off] : a[(uint)(i + 1) * 2 * 4 + off]), 
                                    2) 
                                >> 2);
                }
                for (i = winHX0; i < winHX1; i++)
                {
                    for (off = 0; off < 4; off++)
                        a[(uint)i * 2 * 4 + off] = 
                            MyMath.int_add_no_overflow(
                                a[(uint)i * 2 * 4 + off],
                                MyMath.int_add_no_overflow(
                                    i < 0 ? a[(1 + (uint)0 * 2) * 4 + off] : i >= sn ? a[(1 + (uint)(sn - 1) * 2) * 4 + off] : a[(1 + (uint)i * 2) * 4 + off], 
                                    i - 1 < 0 ? a[(1 + (uint)0 * 2) * 4 + off] : i - 1 >= sn ? a[(1 + (uint)(sn - 1) * 2) * 4 + off] : a[(1 + (uint)(i - 1) * 2) * 4 + off]) 
                                >> 1);
                }
            }
        }
    }

    private static void v8dwt_interleave_h(V4dwtLocal dwt, int[] aAr, int a, int width, uint remainingHeight)
    {
        var biAr = dwt.Wavelet;
        var bi = dwt.Cas * NbEltsV8; //C# NB_ELTS_V8 is because we're indexing into dwt.wavelet
        var x0 = (int)dwt.WinLX0;
        var x1 = (int)dwt.WinLX1;

        //C# impl. note:
        //Workaround for C's ability to treat float and
        //int as raw data.
        var fi = new IntOrFloat();

        for (var k = 0; k < 2; ++k)
        {
            //C# impl note. (a & 0x0f) and (bi & 0x0f) checks if a "pointer" is aligned.
            //Now, on C# the base pointer is always aligned, so this test should still work
            if (remainingHeight >= NbEltsV8 && (a & 0x0f) == 0 && (bi & 0x0f) == 0)
            {
                // Fast code path
                // C# - I've not done any benchmarking.
                for (var i = x0; i < x1; ++i)
                {
                    int j = a + i, dst = bi + i * 2 * NbEltsV8;
                    fi.I = aAr[j];
                    biAr[dst + 0] = fi.F;
                    j += width; fi.I = aAr[j];
                    biAr[dst + 1] = fi.F;
                    j += width; fi.I = aAr[j];
                    biAr[dst + 2] = fi.F;
                    j += width; fi.I = aAr[j];
                    biAr[dst + 3] = fi.F;
                    j += width; fi.I = aAr[j];
                    biAr[dst + 4] = fi.F;
                    j += width; fi.I = aAr[j];
                    biAr[dst + 5] = fi.F;
                    j += width; fi.I = aAr[j];
                    biAr[dst + 6] = fi.F;
                    j += width; fi.I = aAr[j];
                    biAr[dst + 7] = fi.F;
                }
            }
            else
            {
                // Slow code path
                for (var i = x0; i < x1; ++i)
                {
                    int j = a + i, dst = bi + i * 2 * NbEltsV8;
                    fi.I = aAr[j];
                    biAr[dst + 0] = fi.F;
                    j += width;
                    if (remainingHeight == 1) continue;
                    fi.I = aAr[j];
                    biAr[dst + 1] = fi.F;
                    j += width;
                    if (remainingHeight == 2) continue;
                    fi.I = aAr[j];
                    biAr[dst + 2] = fi.F;
                    j += width;
                    if (remainingHeight == 3) continue;
                    fi.I = aAr[j];
                    biAr[dst + 3] = fi.F;
                    j += width;
                    if (remainingHeight == 4) continue;
                    fi.I = aAr[j];
                    biAr[dst + 4] = fi.F;
                    j += width;
                    if (remainingHeight == 5) continue;
                    fi.I = aAr[j];
                    biAr[dst + 5] = fi.F;
                    j += width;
                    if (remainingHeight == 6) continue;
                    fi.I = aAr[j];
                    biAr[dst + 6] = fi.F;
                    j += width;
                    if (remainingHeight == 7) continue;
                    fi.I = aAr[j];
                    biAr[dst + 7] = fi.F;
                }
            }

            bi = (1 - dwt.Cas) * NbEltsV8;
            a += dwt.Sn;
            x0 = (int) dwt.WinHX0;
            x1 = (int) dwt.WinHX1;
        }
    }

    //2.5
    private static void v8dwt_interleave_v(V4dwtLocal dwt, int[] aAr, int a, int width, int nEltsRead)
    {
        //C# Impl. Note that bi_ar is not a float[], but a wavelet array where each entery has 8 floating points.
        //         This while the a array is a plain float[]. I.e how they are to be indexed differers
        var biAr = dwt.Wavelet;
        var bi = dwt.Cas;

        for (var i = (int)dwt.WinLX0; i < dwt.WinLX1; ++i)
            Buffer.BlockCopy(aAr, (a + i * width) * sizeof(float), biAr, (bi + i * 2) * NbEltsV8 * sizeof(float), nEltsRead * sizeof(float));
        a += dwt.Sn * width;
        bi = 1 - dwt.Cas;
        for (var i = (int)dwt.WinHX0; i < dwt.WinHX1; ++i)
            Buffer.BlockCopy(aAr, (a + i * width) * sizeof(float), biAr, (bi + i * 2) * NbEltsV8 * sizeof(float), nEltsRead * sizeof(float));
    }

    //2.5
    private static void v8dwt_interleave_partial_h(V4dwtLocal dwt, SparseArrayInt32 sa, uint saLine, uint remainingHeight)
    {
        for (uint i = 0; i < remainingHeight; i++)
        {
            var ret = sa.read(dwt.WinLX0, saLine + i,
                dwt.WinLX1, saLine + i + 1,
                /* Nasty cast from float* to int32* */
                dwt.Wavelet, (dwt.Cas + 2 * (int)dwt.WinLX0) * NbEltsV8 + (int)i, //C# dwt.wavelet index must be multiplied with NB_ELTS_V8
                2 * NbEltsV8, 0, true);
            Debug.Assert(ret);
            ret = sa.read((uint)dwt.Sn + dwt.WinHX0, saLine + i,
                (uint)dwt.Sn + dwt.WinHX1, saLine + i + 1,
                /* Nasty cast from float* to int32* */
                dwt.Wavelet, (1 - dwt.Cas + 2 * (int)dwt.WinHX0) * NbEltsV8 + (int)i,
                2 * NbEltsV8, 0, true);
            Debug.Assert(ret);
        }
    }

    //2.5
    private static void v8dwt_interleave_partial_v(V4dwtLocal dwt, SparseArrayInt32 sa, uint saCol, uint nbEltsRead)
    {
        var ret = sa.read(saCol, dwt.WinLX0,
            saCol + nbEltsRead, dwt.WinLX1,
            /* Nasty cast from float* to int32* */
            dwt.Wavelet, (dwt.Cas + 2 * (int)dwt.WinLX0) * NbEltsV8,
            1, 2 * NbEltsV8, true);
        Debug.Assert(ret);
        ret = sa.read(saCol, (uint)dwt.Sn + dwt.WinHX0,
            saCol + nbEltsRead, (uint)dwt.Sn + dwt.WinHX1,
            /* Nasty cast from float* to int32* */
            dwt.Wavelet, (1 - dwt.Cas + 2 * (int)dwt.WinHX0) * NbEltsV8,
            1, 2 * NbEltsV8, true);
        Debug.Assert(ret);
    }

    /// <summary>
    /// Forward lazy transform (horizontal).
    /// </summary>
    /// <remarks>2.5 - opj_dwt_deinterleave_h</remarks>
    private static void Deinterleave_h(int[] a, int[] b, int bPt, int dn, int sn, int cas)
    {
        var dest = bPt;
        var src = 0 + cas;

        for (var i = 0; i < sn; i++)
        {
            b[dest++] = a[src];
            src += 2;
        }

        dest = bPt + sn;
        src = 0 + 1 - cas;

        for (var i = 0; i < dn; i++)
        {
            b[dest++] = a[src];
            src += 2;
        }
    }

    /// <summary>
    /// Inverse 5-3 wavelet transform in 2-D
    /// </summary>
    /// <param name="tilec">Tile component information (current tile)</param>
    /// <param name="numres">Number of resolution levels to decode</param>
    /// <remarks>
    /// 2.5 - opj_dwt_decode
    /// </remarks>
    internal static bool Decode(TileCoder tcd, TcdTilecomp tilec, uint numres)
    {
        if (tcd.WholeTileDecoding)
            return decode_tile(tilec, numres, tcd.DisableMultiThreading);
        else
            return DecodePartialTile(tilec, numres);
    }

    /// <summary>
    /// Forward 9-7 wavelet transform in 1-D.
    /// </summary>
    /// <remarks>2.5 - opj_dwt_encode_1_real</remarks>
    private static void Encode_1_real(int[] w, int dn, int sn, int cas)
    {
        int a, b;
        if (cas == 0)
        {
            a = 0;
            b = 1;
        }
        else
        {
            a = 1;
            b = 0;
        }
        EncodeStep2(w, a, b + 1, (uint)dn, (uint)Math.Min(dn, sn - b), DwtAlpha);
        EncodeStep2(w, b, a + 1, (uint)sn, (uint)Math.Min(sn, dn - a), DwtBeta);
        EncodeStep2(w, a, b + 1, (uint)dn, (uint)Math.Min(dn, sn - b), DwtGamma);
        EncodeStep2(w, b, a + 1, (uint)sn, (uint)Math.Min(sn, dn - a), DwtDelta);

        if (a == 0)
        {
            EncodeStep1Combined(w, 0, (uint)sn, (uint)dn, InvK, K);
        }
        else
        {
            EncodeStep1Combined(w, 0, (uint)dn, (uint)sn, K, InvK);
        }
    }

    //2.5 - opj_dwt_encode_and_deinterleave_v_real
    private static void EncodeAndDeinterleaveV_Real(int[] arr, int aPt, int[] tmp, uint height, bool even, uint strideWidth, uint cols)
    {
        if (height == 1)
            return;

        var sn = (int)((height + (even ? 1u : 0u)) >> 1);
        var dn = (int)(height - sn);
        int a, b;

        FetchColsVerticalPass(arr, aPt, tmp, height, (int)strideWidth, cols);

        if (even)
        {
            a = 0;
            b = 1;
        }
        else
        {
            a = 1;
            b = 0;
        }
        v8dwt_encode_step2(tmp, a * NbEltsV8, (b + 1) * NbEltsV8, (uint)dn, (uint)Math.Min(dn, sn - b), DwtAlpha);
        v8dwt_encode_step2(tmp, b * NbEltsV8, (a + 1) * NbEltsV8, (uint)sn, (uint)Math.Min(sn, dn - a), DwtBeta);
        v8dwt_encode_step2(tmp, a * NbEltsV8, (b + 1) * NbEltsV8, (uint)dn, (uint)Math.Min(dn, sn - b), DwtGamma);
        v8dwt_encode_step2(tmp, b * NbEltsV8, (a + 1) * NbEltsV8, (uint)sn, (uint)Math.Min(sn, dn - a), DwtDelta);
        v8dwt_encode_step1(tmp, b * NbEltsV8, (uint)dn, K);
        v8dwt_encode_step1(tmp, a * NbEltsV8, (uint)sn, InvK);

        if (cols == NbEltsV8)
        {
            DeinterleaveV_Cols(tmp, arr, aPt, dn, sn, (int)strideWidth, even ? 0 : 1, NbEltsV8);
        }
        else
        {
            DeinterleaveV_Cols(tmp, arr, aPt, dn, sn, (int)strideWidth, even ? 0 : 1, cols);
        }
    }

    /// <summary>
    /// Forward 5-3 transform, for the vertical pass, processing cols columns
    /// where cols <= NB_ELTS_V8
    /// </summary>
    /// <remarks>2.5 - opj_dwt_encode_and_deinterleave_v</remarks>
    private static void EncodeAndDeinterleaveV(int[] a, int aPt, int[] tmp, uint height, bool even, uint strideWidth, uint cols)
    {
        var sn = (height + (even ? 1u : 0u)) >> 1;
        var dn = height - sn;

        FetchColsVerticalPass(a, aPt, tmp, height, (int)strideWidth, cols);

        //C# Snip SSE2 code

        if (even)
        {
            uint c;
            if (height > 1)
            {
                uint i;
                for (i = 0; i + 1 < sn; i++)
                {
                    for (c = 0; c < 8; c++)
                    {
                        tmp[(1 + i * 2) * 8 + c] -= (tmp[i * 2 * 8 + c] + tmp[(i + 1) * 2 * 8 + c]) >> 1;
                    }
                }
                if (height % 2 == 0)
                {
                    for (c = 0; c < 8; c++)
                    {
                        tmp[(1 + i * 2) * 8 + c] -= tmp[i * 2 * 8 + c];
                    }
                }
                for (c = 0; c < 8; c++)
                {
                    tmp[0 * 2 * 8 + c] += (tmp[(1 + 0 * 2) * 8 + c] + tmp[(1 + 0 * 2) * 8 + c] + 2) >> 2;
                }
                for (i = 1; i < dn; i++)
                {
                    for (c = 0; c < 8; c++)
                    {
                        tmp[i * 2 * 8 + c] += (tmp[(1 + (i - 1) * 2) * 8 + c] + tmp[(1 + i * 2) * 8 + c] + 2) >> 2;
                    }
                }
                if (height % 2 == 1)
                {
                    for (c = 0; c < 8; c++)
                    {
                        tmp[i * 2 * 8 + c] += (tmp[(1 + (i - 1) * 2) * 8 + c] + tmp[(1 + (i - 1) * 2) * 8 + c] + 2) >> 2;
                    }
                }
            }
        }
        else
        {
            uint c;
            if (height == 1)
            {
                for (c = 0; c < 8; c++)
                {
                    tmp[0 * 2 * 8 + c] *= 2;
                }
            }
            else
            {
                uint i;
                for (c = 0; c < 8; c++)
                {
                    tmp[0 * 2 * 8 + c] -= tmp[(1 + 0 * 2) * 8 + c];
                }
                for (i = 1; i < sn; i++)
                {
                    for (c = 0; c < 8; c++)
                    {
                        tmp[i * 2 * 8 + c] -= (tmp[(1 + i * 2) * 8 + c] + tmp[(1 + (i - 1) * 2) * 8 + c]) >> 1;
                    }
                }
                if (height % 2 == 1)
                {
                    for (c = 0; c < 8; c++)
                    {
                        tmp[i * 2 * 8 + c] -= tmp[(1 + (i - 1) * 2) * 8 + c];
                    }
                }
                for (i = 0; i + 1 < dn; i++)
                {
                    for (c = 0; c < 8; c++)
                    {
                        tmp[(1 + i * 2) * 8 + c] += (tmp[i * 2 * 8 + c] + tmp[(i + 1) * 2 * 8 + c] + 2) >> 2;
                    }
                }
                if (height % 2 == 0)
                {
                    for (c = 0; c < 8; c++)
                    {
                        tmp[(1 + i * 2) * 8 + c] += (tmp[i * 2 * 8 + c] + tmp[i * 2 * 8 + c] + 2) >> 2;
                    }
                }
            }
        }

        if (cols == NbEltsV8)
        {
            DeinterleaveV_Cols(tmp, a, aPt, (int)dn, (int)sn,
                (int)strideWidth, even ? 0 : 1, NbEltsV8);
        }
        else
        {
            DeinterleaveV_Cols(tmp, a, aPt, (int)dn, (int)sn,
                (int)strideWidth, even ? 0 : 1, cols);
        }
    }

    /// <summary>
    /// Deinterleave result of forward transform, where cols <= NB_ELTS_V8
    /// and src contains NB_ELTS_V8 consecutive values for up to NB_ELTS_V8
    /// columns
    /// </summary>
    /// <remarks>2.5 - opj_dwt_deinterleave_v_cols</remarks>
    private static void DeinterleaveV_Cols(int[] src, int[] dst, int dstPt, int dn, int sn, int strideWidth, int cas, uint cols)
    {
        var orgDest = dstPt;
        var srcPt = cas * NbEltsV8;
        var i = sn;
        for (var k = 0; k < 2; k++)
        {
            while(i-- != 0)
            {
                if (cols == NbEltsV8)
                {
                    Buffer.BlockCopy(src, srcPt * sizeof(int), dst, dstPt * sizeof(int), NbEltsV8 * sizeof(int));
                }
                else
                {
                    var c = 0;
                    switch(cols)
                    {
                        case 7:
                            dst[dstPt + c] = src[srcPt + c];
                            c++;
                            goto case 6;
                        case 6:
                            dst[dstPt + c] = src[srcPt + c];
                            c++;
                            goto case 5;
                        case 5:
                            dst[dstPt + c] = src[srcPt + c];
                            c++;
                            goto case 4;
                        case 4:
                            dst[dstPt + c] = src[srcPt + c];
                            c++;
                            goto case 3;
                        case 3:
                            dst[dstPt + c] = src[srcPt + c];
                            c++;
                            goto case 2;
                        case 2:
                            dst[dstPt + c] = src[srcPt + c];
                            c++;
                            goto default;
                        default:
                            dst[dstPt + c] = src[srcPt + c];
                            break;
                    }
                }
                dstPt += strideWidth;
                srcPt += 2 * NbEltsV8;
            }

            dstPt = orgDest + sn * strideWidth;
            srcPt = (1 - cas) * NbEltsV8;
            i = dn;
        }
    }

    /// <summary>
    /// Process one line for the horizontal pass of the 9x7 forward transform
    /// </summary>
    /// <remarks>2.5 - opj_dwt_encode_and_deinterleave_h_one_row_real</remarks>
    private static void EncodeAndDeinterleaveH_OneRowReal(int[] row, int rowPt, int[] tmp, uint width, bool even)
    {
        if (width == 1)
            return;
        var sn = (int)((width + (even ? 1 : 0)) >> 1);
        var dn = (int)(width - (uint)sn);
        Buffer.BlockCopy(row, rowPt * sizeof(int), tmp, 0, (int)(width * sizeof(int)));
        Encode_1_real(tmp, dn, sn, even ? 0 : 1);
        Deinterleave_h(tmp, row, rowPt, dn, sn, even ? 0 : 1);
    }

    /// <summary>
    /// Process one line for the horizontal pass of the 5x3 forward transform
    /// </summary>
    /// <remarks>2.5 - opj_dwt_encode_and_deinterleave_h_one_row</remarks>
    private static void EncodeAndDeinterleaveH_OneRow(int[] row, int rowPt, int[] tmp, uint width, bool even)
    {
        var sn = (int)((width + (even ? 1 : 0)) >> 1);
        var dn = (int)(width - (uint)sn);

        if (even)
        {
            if (width > 1)
            {
                int i;
                for (i = 0; i < sn - 1; i++)
                {
                    tmp[sn + i] = row[rowPt + 2 * i + 1] - ((row[rowPt + i * 2] + row[rowPt + (i + 1) * 2]) >> 1);
                }
                if (width % 2 == 0)
                {
                    tmp[sn + i] = row[rowPt + 2 * i + 1] - row[rowPt + i * 2];
                }
                row[rowPt + 0] += (tmp[sn] + tmp[sn] + 2) >> 2;
                for (i = 1; i < dn; i++)
                {
                    row[rowPt + i] = row[rowPt + 2 * i] + ((tmp[sn + (i - 1)] + tmp[sn + i] + 2) >> 2);
                }
                if (width % 2 == 1)
                {
                    row[rowPt + i] = row[rowPt + 2 * i] + ((tmp[sn + (i - 1)] + tmp[sn + (i - 1)] + 2) >> 2);
                }
                Buffer.BlockCopy(tmp, sn * sizeof(int), row, (rowPt + sn) * sizeof(int), dn * sizeof(int));
            }
        }
        else
        {
            if (width == 1)
            {
                row[rowPt] *= 2;
            }
            else
            {
                int i;
                tmp[sn + 0] = row[rowPt + 0] - row[rowPt + 1];
                for (i = 1; i < sn; i++)
                {
                    tmp[sn + i] = row[rowPt + 2 * i] - ((row[rowPt + 2 * i + 1] + row[rowPt + 2 * (i - 1) + 1]) >> 1);
                }
                if (width % 2 == 1)
                {
                    tmp[sn + i] = row[rowPt + 2 * i] - row[rowPt + 2 * (i - 1) + 1];
                }

                for (i = 0; i < dn - 1; i++)
                {
                    row[rowPt + i] = row[rowPt + 2 * i + 1] + ((tmp[sn + i] + tmp[sn + i + 1] + 2) >> 2);
                }
                if (width % 2 == 0)
                {
                    row[rowPt + i] = row[rowPt + 2 * i + 1] + ((tmp[sn + i] + tmp[sn + i] + 2) >> 2);
                }
                Buffer.BlockCopy(tmp, sn * sizeof(int), row, (rowPt + sn) * sizeof(int), dn * sizeof(int));
            }
        }
    }

    /** Fetch up to cols <= NB_ELTS_V8 for each line, and put them in tmpOut */
    /* that has a NB_ELTS_V8 interleave factor. */
    //2.5
    private static void FetchColsVerticalPass(int[] a, int array, int[] tmp, uint height, int strideWidth, uint cols)
    {
        if (cols == NbEltsV8) {
            for (var k = 0; k < height; ++k) {
                Buffer.BlockCopy(a, (array + k * strideWidth) * sizeof(int), tmp, NbEltsV8 * k * sizeof(int), NbEltsV8 * sizeof(int));
            }
        } else {
            for (var k = 0; k < height; ++k) {
                var c = 0;
                for (; c < cols; c++) {
                    tmp[NbEltsV8 * k + c] = a[array + c + k * strideWidth];
                }
                for (; c < NbEltsV8; c++) {
                    tmp[NbEltsV8 * k + c] = 0;
                }
            }
        }
    }

    //2.5
    private static void decode_h_func(DecodeHJob job)
    {
        for(var j = (int)job.MinJ; j < job.MaxJ; j++)
        {
            idwt53_h(job.H, job.Tiled, job.Tiledp + j * (int)job.W);
        }
    }

    //2.5
    private static void dwt97_decode_h_func(Dwt97DecodeHJob job)
    {
        var w = (int)job.W;
        var ajAr = job.Aj;
        var aj = job.Ajp;
        var fi = new IntOrFloat();

        for (var j = 0; j + NbEltsV8 <= job.NbRows; j += NbEltsV8)
        {
            v8dwt_interleave_h(job.H, ajAr, aj, w, NbEltsV8);
            v8dwt_decode(job.H);

            // To be adapted if NB_ELTS_V8 changes
            for (var k = 0; k < job.Rw; k++)
            {
                //C# note: Org. impl stores the wavlet as a struct with four
                //floating points. Here it's stored as a continious array. 
                var kWavelet = k * NbEltsV8;

                fi.F = job.H.Wavelet[kWavelet + 0];
                ajAr[aj + k] = fi.I;
                fi.F = job.H.Wavelet[kWavelet + 1];
                ajAr[aj + k + w] = fi.I;
                fi.F = job.H.Wavelet[kWavelet + 2];
                ajAr[aj + k + w * 2] = fi.I;
                fi.F = job.H.Wavelet[kWavelet + 3];
                ajAr[aj + k + w * 3] = fi.I;
            }
            for (var k = 0; k < job.Rw; k++)
            {
                var kWavelet = k * NbEltsV8;

                fi.F = job.H.Wavelet[kWavelet + 4];
                ajAr[aj + k + w * 4] = fi.I;
                fi.F = job.H.Wavelet[kWavelet + 5];
                ajAr[aj + k + w * 5] = fi.I;
                fi.F = job.H.Wavelet[kWavelet + 6];
                ajAr[aj + k + w * 6] = fi.I;
                fi.F = job.H.Wavelet[kWavelet + 7];
                ajAr[aj + k + w * 7] = fi.I;
            }

            aj += w * NbEltsV8;
        }
    }

    //2.5
    private static void dwt97_decode_v_func(Dwt97DecodeVJob job)
    {
        var ajAr = job.Aj;
        var aj = job.Ajp;

        for (uint j = 0; j + NbEltsV8 <= job.NbColumns; j += NbEltsV8)
        {
            v8dwt_interleave_v(job.V, ajAr, aj, (int)job.W, NbEltsV8);
            v8dwt_decode(job.V);
            for (var k = 0; k < job.Rh; ++k)
                Buffer.BlockCopy(job.V.Wavelet, k * NbEltsV8 * sizeof(float), ajAr, (aj + k * (int)job.W) * sizeof(float), NbEltsV8 * sizeof(float));

            aj += NbEltsV8;
        }
    }

    /// <summary>
    /// Inverse 5-3 wavelet transform in 1-D for one row.
    /// </summary>
    /// <remarks>
    /// 2.5
    /// Performs interleave, inverse wavelet transform and copy back to buffer
    /// </remarks>
    private static void idwt53_h(DwtLocal dwt, int[] tiled, int tiledp)
    {
#if STANDARD_SLOW_VERSION
            Interleave_h(dwt, tiledp, tiled);
            Decode_1(dwt);
            Buffer.BlockCopy(dwt.mem, 0, tiled, tiledp * sizeof(int), (dwt.sn + dwt.dn) * sizeof(int));
#else
        var sn = dwt.Sn;
        var len = sn + dwt.Dn;
        if (dwt.Cas == 0)
        { /* Left-most sample is on even coordinate */
            if (len > 1)
            {
                idwt53_h_cas0(dwt.Mem, sn, len, tiled, tiledp);
            }
            else
            {
                /* Unmodified value */
            }
        }
        else
        { /* Left-most sample is on odd coordinate */
            if (len == 1)
            {
                tiled[tiledp] /= 2;
            }
            else if (len == 2)
            {
                var o = dwt.Mem;
                var inEven = tiledp + sn;
                var inOdd = tiledp;
                o[1] = tiled[inOdd] - ((tiled[inEven] + 1) >> 1);
                o[0] = tiled[inEven] + o[1];
                Buffer.BlockCopy(dwt.Mem, 0, tiled, tiledp * sizeof(int), len * sizeof(int));
            }
            else if (len > 2)
            {
                opj_idwt53_h_cas1(dwt.Mem, sn, len, tiled, tiledp);
            }
        }
#endif
    }

#if STANDARD_SLOW_VERSION
        /// <summary>
        /// Inverse 5-3 wavelet transform in 1-D.
        /// </summary>
        //2.1
        static void Decode_1(dwt_local v)
        {
            Decode_1_(v.mem, v.dn, v.sn, v.cas);
        }

        /// <summary>
        /// Inverse 5-3 wavelet transform in 1-D.
        /// </summary>
        //2.1
        static void Decode_1_(int[] a, int dn, int sn, int cas)
        {
            if (cas == 0)
            {
                if (dn > 0 || sn > 1)
                {
                    for (int i = 0; i < sn; i++)
                    {
                        //C# impl. of the macro: D_(i - 1) + D_(i)
                        int D1 = (i - 1) < 0 ? 0 : (i - 1) >= dn ? (dn - 1) : (i - 1);
                        int D2 = i < 0 ? 0 : i >= dn ? (dn - 1) : i;

                        //S(i) -= (D_(i - 1) + D_(i) + 2) >> 2;
                        a[i * 2] -= (a[1 + D1 * 2] + a[1 + D2 * 2] + 2) >> 2;
                    }
                    for (int i = 0; i < dn; i++)
                    {
                        //C# impl. of the macro: S_(i) + S_(i + 1)
                        int S1 = i < 0 ? 0 : i >= sn ? (sn - 1) : i;
                        int S2 = (i + 1) < 0 ? 0 : (i + 1) >= sn ? (sn - 1) : (i + 1);

                        //D(i) += (S_(i) + S_(i + 1)) >> 1;
                        a[1 + i * 2] += (a[S1 * 2] + a[S2 * 2]) >> 1;
                    }
                    //#define S(i) a[(i)*2]
                    //#define D(i) a[(1+(i)*2)]
                     //#define S_(i) ((i)<0?S(0):((i)>=sn?S(sn-1):S(i)))
                    //#define SS_(i) ((i)<0?S(0):((i)>=dn?S(dn-1):S(i)))
                     //#define D_(i) ((i)<0?D(0):((i)>=dn?D(dn-1):D(i)))
                    //#define DD_(i) ((i)<0?D(0):((i)>=sn?D(sn-1):D(i)))
                }
            }
            else
            {
                if (sn == 0 && dn == 1)
                    a[0] /= 2;
                else
                {
                    for (int i = 0; i < sn; i++)
                    {
                        //C# impl. of the macro: SS_(i) + SS_(i + 1)
                        int SS1 = i < 0 ? 0 : i >= dn ? (dn - 1) : i;
                        int SS2 = (i + 1) < 0 ? 0 : (i + 1) >= dn ? (dn - 1) : (i + 1);

                        //D(i) -= (SS_(i) + SS_(i + 1) + 2) >> 2;
                        a[1 + i * 2] -= (a[SS1 * 2] + a[SS2 * 2] + 2) >> 2;
                    }
                    for (int i = 0; i < dn; i++)
                    {
                        //C# impl. of the macro: DD_(i) + DD_(i - 1)
                        int DD1 = i < 0 ? 0 : i >= sn ? (sn - 1) : i;
                        int DD2 = (i - 1) < 0 ? 0 : (i - 1) >= sn ? (sn - 1) : (i - 1);

                        //S(i) += (DD_(i) + DD_(i - 1)) >> 1;
                        a[i * 2] += (a[1 + DD1 * 2] + a[1 + DD2 * 2]) >> 1;
                    }
                }
            }
        }

        
        /// <summary>
        /// Inverse lazy transform (vertical).
        /// </summary>
        //2.1
        static void Interleave_v(dwt_local v, int a, int[] a_ar, int x)
        {
            int ai = a;
            int bi = v.cas;
            int[] b_ar = v.mem;
            int i = v.sn;
            while (i-- != 0)
            {
                b_ar[bi] = a_ar[ai];
                bi += 2;
                ai += x;
            }
            ai = a + (v.sn * x);
            bi = 1 - v.cas;
            i = v.dn;
            while (i-- != 0)
            {
                b_ar[bi] = a_ar[ai];
                bi += 2;
                ai += x;
            }
        }
#else
    //2.5
    private static void idwt53_h_cas0(int[] tmp, int sn, int len, int[] tiled, int tiledp)
    {
        Debug.Assert(len > 1);

        int i, j;
        var inEven = tiledp;
        var inOdd = tiledp + sn;

#if TWO_PASS_VERSION
            /* For documentation purpose: performs lifting in two iterations, */
            /* but without explicit interleaving */

            /* Even */
            tmp[0] = tiled[in_even] - ((tiled[in_odd] + 1) >> 1);
            for (i = 2, j = 0; i <= len - 2; i += 2, j++)
            {
                tmp[i] = tiled[in_even + j + 1] - ((tiled[in_odd + j] + tiled[in_odd + j + 1] + 2) >> 2);
            }
            if ((len & 1) != 0)
            { /* if len is odd */
                tmp[len - 1] = tiled[in_even + (len - 1) / 2] - ((tiled[in_odd + (len - 2) / 2] + 1) >> 1);
            }

            /* Odd */
            for (i = 1, j = 0; i < len - 1; i += 2, j++)
            {
                tmp[i] = tiled[in_odd + j] + ((tmp[i - 1] + tmp[i + 1]) >> 1);
            }
            if ((len & 1) == 0)
            { /* if len is even */
                tmp[len - 1] = tiled[in_odd + (len - 1) / 2] + tmp[len - 2];
            }
#else
        // Improved version of the TWO_PASS_VERSION:
        // Performs lifting in one single iteration. Saves memory
        // accesses and explicit interleaving.

        var s1N = tiled[inEven];
        var d1N = tiled[inOdd];
        var s0N = s1N - ((d1N + 1) >> 1);

        for (i = 0, j = 1; i < len - 3; i += 2, j++)
        {
            var d1C = d1N;
            var s0C = s0N;

            s1N = tiled[inEven + j];
            d1N = tiled[inOdd + j];

            s0N = s1N - ((d1C + d1N + 2) >> 2);

            tmp[i] = s0C;
            tmp[i + 1] = MyMath.int_add_no_overflow(d1C, MyMath.int_add_no_overflow(s0C,
                s0N) >> 1);
        }

        tmp[i] = s0N;

        if ((len & 1) != 0)
        {
            tmp[len - 1] = tiled[inEven + (len - 1) / 2] - ((d1N + 1) >> 1);
            tmp[len - 2] = d1N + ((s0N + tmp[len - 1]) >> 1);
        }
        else
        {
            tmp[len - 1] = d1N + s0N;
        }
#endif
        Buffer.BlockCopy(tmp, 0, tiled, tiledp * sizeof(int), len * sizeof(int));
    }

    //2.5
    private static void opj_idwt53_h_cas1(int[] tmp, int sn, int len, int[] tiled, int tiledp)
    {
        Debug.Assert(len > 2);

        int i, j;
        var inEven = tiledp + sn;
        var inOdd = tiledp;

#if TWO_PASS_VERSION
            /* For documentation purpose: performs lifting in two iterations, */
            /* but without explicit interleaving */

            /* Odd */
            for (i = 1, j = 0; i < len - 1; i += 2, j++)
            {
                tmp[i] = tiled[in_odd + j] - ((tiled[in_even + j] + tiled[in_even + j + 1] + 2) >> 2);
            }
            if ((len & 1) == 0)
            {
                tmp[len - 1] = tiled[in_odd + len / 2 - 1] - ((tiled[in_even + len / 2 - 1] + 1) >> 1);
            }

            /* Even */
            tmp[0] = tiled[in_even] + tmp[1];
            for (i = 2, j = 1; i < len - 1; i += 2, j++)
            {
                tmp[i] = tiled[in_even + j] + ((tmp[i + 1] + tmp[i - 1]) >> 1);
            }
            if ((len & 1) != 0)
            {
                tmp[len - 1] = tiled[in_even + len / 2] + tmp[len - 2];
            }
#else
        int dn;

        /* Improved version of the TWO_PASS_VERSION: */
        /* Performs lifting in one single iteration. Saves memory */
        /* accesses and explicit interleaving. */
        var s1 = tiled[inEven + 1];
        var dc = tiled[inOdd] - ((tiled[inEven] + s1 + 2) >> 2);
        tmp[0] = tiled[inEven] + dc;

        var end = len - 2 - ((len & 1) == 0 ? 1 : 0);
        for (i = 1, j = 1; i < end; i += 2, j++)
        {

            var s2 = tiled[inEven + j + 1];

            dn = tiled[inOdd + j] - ((s1 + s2 + 2) >> 2);
            tmp[i] = dc;
            tmp[i + 1] = MyMath.int_add_no_overflow(s1, MyMath.int_add_no_overflow(dn, dc) >> 1);

            dc = dn;
            s1 = s2;
        }

        tmp[i] = dc;

        if ((len & 1) == 0)
        {
            dn = tiled[inOdd + len / 2 - 1] - ((s1 + 1) >> 1);
            tmp[len - 2] = s1 + ((dn + dc) >> 1);
            tmp[len - 1] = dn;
        }
        else
        {
            tmp[len - 1] = s1 + dc;
        }
#endif
        Buffer.BlockCopy(tmp, 0, tiled, tiledp * sizeof(int), len * sizeof(int));
    }

    //2.5
    private static void idwt3_v_cas0(int[] tmp, int sn, int len, int[] tiled, int tiledpCol, int stride)
    {
        int i, j;

        Debug.Assert(len > 1);

        /* Performs lifting in one single iteration. Saves memory */
        /* accesses and explicit interleaving. */
        var s1N = tiled[tiledpCol];
        var d1N = tiled[tiledpCol + sn * stride];
        var s0N = s1N - ((d1N + 1) >> 1);

        for (i = 0, j = 0; i < len - 3; i += 2, j++)
        {
            var d1C = d1N;
            var s0C = s0N;

            s1N = tiled[tiledpCol + (j + 1) * stride];
            d1N = tiled[tiledpCol + (sn + j + 1) * stride];

            s0N = MyMath.int_sub_no_overflow(s1N,
                MyMath.int_add_no_overflow(MyMath.int_add_no_overflow(d1C, d1N), 2) >> 2);

            tmp[i] = s0C;
            tmp[i + 1] = MyMath.int_add_no_overflow(d1C, MyMath.int_add_no_overflow(s0C,
                s0N) >> 1);
        }

        tmp[i] = s0N;

        if ((len & 1) != 0)
        {
            tmp[len - 1] =
                tiled[tiledpCol + (len - 1) / 2 * stride] -
                ((d1N + 1) >> 1);
            tmp[len - 2] = d1N + ((s0N + tmp[len - 1]) >> 1);
        }
        else
        {
            tmp[len - 1] = d1N + s0N;
        }

        for (i = 0; i < len; ++i)
        {
            tiled[tiledpCol + i * stride] = tmp[i];
            //if (214928 == tiledp_col + i * stride)
            //{
            //    Debug.Write("Hello");
            //}
        }
    }

    //2.5
    private static void idwt3_v_cas1(int[] tmp, int sn, int len, int[] tiled, int tiledpCol, int stride)
    {
        int i, j;
        int dn;
        var inEven = tiledpCol + sn * stride;
        var inOdd = tiledpCol;

        Debug.Assert(len > 2);

        // Performs lifting in one single iteration. Saves memory
        // accesses and explicit interleaving.
        var s1 = tiled[inEven + stride];
        var dc = tiled[inOdd] - ((tiled[inEven] + s1 + 2) >> 2);
        tmp[0] = tiled[inEven] + dc;
        var end = len - 2 - ((len & 1) == 0 ? 1 : 0);
        for (i = 1, j = 1; i < end; i += 2, j++)
        {

            var s2 = tiled[inEven + (j + 1) * stride];

            dn = tiled[inOdd + j * stride] - ((s1 + s2 + 2) >> 2);
            tmp[i] = dc;
            tmp[i + 1] = s1 + ((dn + dc) >> 1);

            dc = dn;
            s1 = s2;
        }
        tmp[i] = dc;
        if ((len & 1) == 0)
        {
            dn = tiled[inOdd + (len / 2 - 1) * stride] - ((s1 + 1) >> 1);
            tmp[len - 2] = s1 + ((dn + dc) >> 1);
            tmp[len - 1] = dn;
        }
        else
        {
            tmp[len - 1] = s1 + dc;
        }

        for (i = 0; i < len; ++i)
        {
            tiled[tiledpCol + i * stride] = tmp[i];
        }
    }
#endif

    //2.5
    private static void decode_v_func(DecodeVJob job)
    {
        int j;
        for (j = (int)job.MinJ; j + ParallelCols53 <= job.MaxJ; j += (int)ParallelCols53)
        {
            idwt53_v(job.V, job.Tiled, job.Tiledp +j, (int)job.W, (int)ParallelCols53);
        }
        if (j < job.MaxJ)
            idwt53_v(job.V, job.Tiled, job.Tiledp + j, (int)job.W, (int)(job.MaxJ - j));
    }

    //2.5
    private static void idwt53_v(DwtLocal dwt, int[] tiled, int tiledpCol, int stride, int nbCols)
    {
#if STANDARD_SLOW_VERSION
            for (int c = 0; c < nb_cols; c ++) {
                Interleave_v(dwt, tiledp_col + c, tiled, stride);
                Decode_1(dwt);
                for (int k = 0; k < dwt.sn + dwt.dn; ++k) {
                    tiled[tiledp_col + c + k * stride] = dwt.mem[k];
                }
            }
#else
        //var mem_copy = new int[dwt.mem.Length];
        //Array.Copy(dwt.mem, mem_copy, dwt.mem.Length);
        //var tiled_copy = new int[tiled.Length];
        //Array.Copy(tiled, tiled_copy, tiled_copy.Length);

        //for (int c = 0; c < nb_cols; c++)
        //{
        //    Interleave_v(dwt, tiledp_col + c, tiled, stride);
        //    Decode_1(dwt);
        //    for (int k = 0; k < dwt.sn + dwt.dn; ++k)
        //    {
        //        tiled[tiledp_col + c + k * stride] = dwt.mem[k];
        //    }
        //}

        //var correct_tiled = tiled;
        //tiled = tiled_copy;
        //dwt.mem = mem_copy;

        var sn = dwt.Sn;
        var len = sn + dwt.Dn;
        if (dwt.Cas == 0)
        {
            //C# Snip SSE2

            if (len > 1)
            {
                for (var c = 0; c < nbCols; c++, tiledpCol++)
                {
                    idwt3_v_cas0(dwt.Mem, sn, len, tiled, tiledpCol, stride);
                }

                ////Check for correctnes
                //for (int c = 0; c < correct_tiled.Length; c++)
                //{
                //    if (correct_tiled[c] != tiled[c])
                //    {
                //        Debug.Write("Nah");
                //    }
                //}
                return;
            }
        }
        else
        {
            if (len == 1)
            {
                for (var c = 0; c < nbCols; c++, tiledpCol++)
                {
                    tiled[tiledpCol] /= 2;
                }
                return;
            }

            if (len == 2)
            {
                var o = dwt.Mem;
                for (var c = 0; c < nbCols; c++, tiledpCol++)
                {
                    var inEven = tiledpCol + sn * stride;
                    var inOdd = tiledpCol;

                    o[1] = tiled[inOdd] - ((tiled[inEven] + 1) >> 1);
                    o[0] = tiled[inEven] + o[1];

                    for (var i = 0; i < len; ++i)
                    {
                        tiled[tiledpCol + i * stride] = o[i] ;
                    }
                }

                return;
            }

            //C# Snip SSE2

            if (len > 2)
            {
                for (var c = 0; c < nbCols; c++, tiledpCol++)
                {
                    idwt3_v_cas1(dwt.Mem, sn, len, tiled, tiledpCol, stride);
                }
                return;
            }
        }
#endif
    }

    /// <summary>
    /// Get norm of 9-7 wavelet.
    /// </summary>
    /// <remarks>2.5 - opj_dwt_getnorm_real</remarks>
    internal static double Getnorm_real(uint level, uint orient)
    {
        if (orient == 0 && level >= 10)
            level = 9;
        else if (orient > 0 && level >= 9)
            level = 8;
        return DwtNormsReal[orient][level];
    }

    /// <summary>
    /// Get norm of 5-3 wavelet.
    /// </summary>
    /// <remarks>2.5 - opj_dwt_getnorm</remarks>
    internal static double Getnorm(uint level, uint orient)
    {
        //FIXME ! This is just a band-aid to avoid a buffer overflow
        if (orient == 0 && level >= 10)
            level = 9;
        else if (orient > 0 && level >= 9)
            level = 8;

        return DwtNorms[orient][level];
    }

    internal delegate int GetGainFunc(int orient);

    //2.5
    private static void GetBandCoordinates(TcdTilecomp tilec,
        uint resno,
        uint bandno,
        uint tcx0,
        uint tcy0,
        uint tcx1,
        uint tcy1,
        out uint tbx0,
        out uint tby0,
        out uint tbx1,
        out uint tby1)
    {
        // Compute number of decomposition for this band. See table F-1
        var nb = resno == 0 ?
            (int)tilec.numresolutions - 1 :
            (int)(tilec.numresolutions - resno);
        /* Map above tile-based coordinates to sub-band-based coordinates per */
        /* equation B-15 of the standard */
        var x0B = bandno & 1;
        var y0B = bandno >> 1;
        //if (tbx0)
        {
            tbx0 = nb == 0 ? tcx0 :
                tcx0 <= (1U << nb - 1) * x0B ? 0 :
                MyMath.uint_ceildivpow2(tcx0 - (1U << (nb - 1)) * x0B, nb);
        }
        //if (tby0)
        {
            tby0 = nb == 0 ? tcy0 :
                tcy0 <= (1U << (nb - 1)) * y0B ? 0 :
                MyMath.uint_ceildivpow2(tcy0 - (1U << (nb - 1)) * y0B, nb);
        }
        //if (tbx1)
        {
            tbx1 = nb == 0 ? tcx1 :
                tcx1 <= (1U << (nb - 1)) * x0B ? 0 :
                MyMath.uint_ceildivpow2(tcx1 - (1U << (nb - 1)) * x0B, nb);
        }
        //if (tby1)
        {
            tby1 = nb == 0 ? tcy1 :
                tcy1 <= (1U << (nb - 1)) * y0B ? 0 :
                MyMath.uint_ceildivpow2(tcy1 - (1U << (nb - 1)) * y0B, nb);
        }
    }

    //2.5 - opj_dwt_segment_grow
    private static void SegmentGrow(uint filterWidth, uint maxSize, ref uint start, ref uint end)
    {
        start = MyMath.uint_subs(start, filterWidth);
        end = MyMath.uint_adds(end, filterWidth);
        end = Math.Min(end, maxSize);
    }

    private delegate void DwgAction(DwtLocal dwt);

    private class DwtLocal
    {
        internal int[] Mem;

        /// <summary>
        /// Number of elements in high pass band
        /// </summary>
        internal int Dn;

        /// <summary>
        /// Number of elements in low pass band
        /// </summary>
        internal int Sn;

        /// <summary>
        /// 0 = start on even coord, 1 = start on odd coord
        /// </summary>
        internal int Cas;

        public DwtLocal Clone()
        {
            return (DwtLocal)MemberwiseClone();
        }
    }

    private class V4dwtLocal
    {
        /// <summary>
        /// Each wavelet is 4 floating points.
        /// </summary>
        /// <remarks>
        ///C# note: Org. impl stores the wavlet as a struct with eight
        ///floating points. Here it's stored as a continious array. 
        ///
        /// This so to make it possible to use Buffer.BlockCopy to
        /// copy the values. 
        ///</remarks>
        internal float[] Wavelet;

        /// <summary>
        /// Number of elements in high pass band
        /// </summary>
        internal int Dn;

        /// <summary>
        /// Number of elements in low pass band
        /// </summary>
        internal int Sn;

        /// <summary>
        ///  0 = start on even coord, 1 = start on odd coord
        /// </summary>
        internal int Cas;

        /// <summary>
        /// Start coord in low pass band
        /// </summary>
        internal uint WinLX0;

        /// <summary>
        /// End coord in low pass band
        /// </summary>
        internal uint WinLX1;

        /// <summary>
        /// Start coord in high pass band
        /// </summary>
        internal uint WinHX0;

        /// <summary>
        /// End coord in high pass band
        /// </summary>
        internal uint WinHX1;

        public V4dwtLocal Clone() { return (V4dwtLocal)MemberwiseClone(); }
    }

    /// <summary>
    /// Used by the MT impl
    /// </summary>
    private class DecodeHJob
    {
        public readonly DwtLocal H;
        public readonly uint Rw;
        public readonly uint W;
        public readonly int[] Tiled;
        public readonly int Tiledp;
        public readonly uint MinJ;
        public readonly uint MaxJ;

        public DecodeHJob(
            DwtLocal h,
            uint rw, uint w,
            int[] tiled, int tiledp,
            uint minJ, uint maxJ
        )
        {
            this.H = h;
            this.Rw = rw;
            this.W = w;
            this.Tiled = tiled;
            this.Tiledp = tiledp;
            this.MinJ = minJ;
            this.MaxJ = maxJ;
        }
    }

    /// <summary>
    /// Used by the MT impl
    /// </summary>
    private class EncodeVJob
    {
        public readonly DwtLocal V;
        public readonly uint Rh;
        public readonly uint W;
        public readonly int[] Tiled;
        public readonly int Tiledp;
        public readonly uint MinJ;
        public readonly uint MaxJ;
        public readonly EncodeAndDeinterleaveVfunc EncodeAndDeinterleaveV;

        public EncodeVJob(
            DwtLocal v,
            uint rh, uint w,
            int[] tiled, int tiledp,
            uint minJ, uint maxJ,
            EncodeAndDeinterleaveVfunc encodeAndDeinterleaveV
        )
        {
            this.V = v;
            this.Rh = rh;
            this.W = w;
            this.Tiled = tiled;
            this.Tiledp = tiledp;
            this.MinJ = minJ;
            this.MaxJ = maxJ;
            this.EncodeAndDeinterleaveV = encodeAndDeinterleaveV;
        }
    }

    /// <summary>
    /// Used by the MT impl
    /// </summary>
    private class EncodeHJob
    {
        public readonly DwtLocal H;
        public readonly uint Rw;
        public readonly uint W;
        public readonly int[] Tiled;
        public readonly int Tiledp;
        public readonly uint MinJ;
        public readonly uint MaxJ;
        public readonly EncodeAndDeinterleaveHOneRowfunc Fn;

        public EncodeHJob(
            DwtLocal h,
            uint rw, uint w,
            int[] tiled, int tiledp,
            uint minJ, uint maxJ,
            EncodeAndDeinterleaveHOneRowfunc fn
        )
        {
            this.H = h;
            this.Rw = rw;
            this.W = w;
            this.Tiled = tiled;
            this.Tiledp = tiledp;
            this.MinJ = minJ;
            this.MaxJ = maxJ;
            this.Fn = fn;
        }
    }

    /// <summary>
    /// Used by the MT impl
    /// </summary>
    private class DecodeVJob
    {
        public readonly DwtLocal V;
        public readonly uint Rh;
        public readonly uint W;
        public readonly int[] Tiled;
        public readonly int Tiledp;
        public readonly uint MinJ;
        public readonly uint MaxJ;

        public DecodeVJob(
            DwtLocal v,
            uint rh, uint w,
            int[] tiled, int tiledp,
            uint minJ, uint maxJ
        )
        {
            this.V = v;
            this.Rh = rh;
            this.W = w;
            this.Tiled = tiled;
            this.Tiledp = tiledp;
            this.MinJ = minJ;
            this.MaxJ = maxJ;
        }
    }

    /// <summary>
    /// Used by the MT impl
    /// </summary>
    private class Dwt97DecodeHJob
    {
        public readonly V4dwtLocal H;
        public readonly uint Rw;
        public readonly uint W;
        public readonly int[] Aj;
        public readonly int Ajp;
        public readonly uint NbRows;

        public Dwt97DecodeHJob(
            V4dwtLocal h,
            uint rw, uint w,
            int[] aj, int ajp,
            uint nbRows
        )
        {
            this.H = h;
            this.Rw = rw;
            this.W = w;
            this.Aj = aj;
            this.Ajp = ajp;
            this.NbRows = nbRows;
        }
    }

    /// <summary>
    /// Used by the MT impl
    /// </summary>
    private class Dwt97DecodeVJob
    {
        public readonly V4dwtLocal V;
        public readonly uint Rh;
        public readonly uint W;
        public readonly int[] Aj;
        public readonly int Ajp;
        public readonly uint NbColumns;

        public Dwt97DecodeVJob(
            V4dwtLocal v,
            uint rh, uint w,
            int[] aj, int ajp,
            uint nbColumns
        )
        {
            this.V = v;
            this.Rh = rh;
            this.W = w;
            this.Aj = aj;
            this.Ajp = ajp;
            this.NbColumns = nbColumns;
        }
    }
}