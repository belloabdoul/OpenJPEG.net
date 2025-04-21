#region License
/*
 * Copyright (c) 2002-2007, Communications and Remote Sensing Laboratory, Universite catholique de Louvain (UCL), Belgium
 * Copyright (c) 2002-2007, Professor Benoit Macq
 * Copyright (c) 2001-2003, David Janssens
 * Copyright (c) 2002-2003, Yannick Verschueren
 * Copyright (c) 2003-2007, Francois-Olivier Devaux and Antonin Descampe
 * Copyright (c) 2005, Herve Drolon, FreeImage Team
 * Copyright (c) 2006-2007, Parvatha Elangovan
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

using System;
using System.IO;
using System.Diagnostics;
using OpenJpeg.Internal;

namespace OpenJpeg;

/// <summary>
/// JPEG-2000 file format reader/writer
/// </summary>
/// <remarks>
/// V.2.1 does things a bit differently. Basically,
/// methods that are to be executed are put into a
/// _procedure_list, then that list is executed.
/// 
/// I don't see the point. Sticking to the old v.1.4
/// way of doing this.
/// </remarks>
internal sealed class Jp2
{
    #region Variables and properties

    /// <summary>
    /// The parent compression info obj.
    /// </summary>
    private readonly CompressionInfo _cinfo;

    /// <summary>
    /// Code stream codex
    /// </summary>
    private readonly J2K _j2K;

    //List<ProcedureDlg> _validation_list;
    //List<ProcedureDlg> _procedure_list;

    private JP2_STATE _state;
    private JP2_IMG_STATE _imgState;

    /// <summary>
    /// Width of the image
    /// </summary>
    private uint _w;

    /// <summary>
    /// Height of the image
    /// </summary>
    private uint _h;

    /// <summary>
    /// Number of componets in the image
    /// </summary>
    private uint _numcomps;

    /// <summary>
    /// Bits per component
    /// </summary>
    private uint _bpc;

    /// <summary>
    /// ColorSpecMethod
    /// </summary>
    private uint _meth;

    /// <summary>
    /// ColorSpace
    /// </summary>
    private uint _enumcs;

    private uint _c, _approx, _precedence, _minversion;
    private JP2_Marker _brand;

    /// <summary>
    /// Unknown color space
    /// </summary>
    private bool _unkC;

    /// <summary>
    /// Intellectual Property
    /// </summary>
    private bool _ipr;

    private JP2_Marker[] _cl;

    private long _j2KCodestreamOffset;

    private JP2Comps[] _comps;

    private Cio _cio;

    /// <summary>
    /// If the image being decoded has a ICC profile, it will be temporarily stored
    /// here. 
    /// </summary>
    private JP2Color _color;

    private bool _hasIhdr;
    private bool HasJp2H => (_state & JP2_STATE.HEADER) != 0;

    #endregion

    #region Init

    internal Jp2(CompressionInfo cinfo, J2K j2K)
    {
        _cinfo = cinfo;
        _j2K = j2K;
            
        //_validation_list = new List<ProcedureDlg>(1);
        //_procedure_list = new List<ProcedureDlg>();
    }

    //2.5.3 - opj_jp2_setup_encoder
    internal bool SetupEncoder(CompressionParameters parameters, JPXImage image)
    {
        if (parameters == null || image == null)
            return false;

        //
        // Sets up the J2K codec
        //

        // Checks if the number of components respects the standard
        if (image.numcomps < 1 || image.numcomps > 16384)
        {
            _cinfo.Error("Invalid number of components specified while setting up JP2 encoder");
            return false;
        }

        _j2K.SetupEncoder(parameters, image);

        //
        // Sets up the JP2 codec
        //

        // Profile box
        _brand = JP2_Marker.JP2;
        _minversion = 0;
        //numcl = _cl.Length
        _cl = new JP2_Marker[1];
        _cl[0] = _brand;

        // Image Header box 
        _numcomps = image.numcomps;
        _comps = GC.AllocateUninitializedArray<JP2Comps>((int)_numcomps);
        _h = image.y1 - image.y0;
        _w = image.x1 - image.x0;
        //Setting bits per componet
        var depth0 = image.comps[0].prec - 1;
        var sign = image.comps[0].sgnd ? 1u : 0u;
        _bpc = depth0 + (sign << 7);
        for (var i = 0; i < image.numcomps; i++)
        {
            var depth = image.comps[i].prec - 1;
            sign = image.comps[i].sgnd ? 1u : 0u;

            //If the pits per component aren't uniform,
            //bpc is set to 255 to signal that.
            if (depth0 != depth)
                _bpc = 255;
        }
        _c = 7;
        _unkC = false;
        _ipr = false;

        //BitsPerComponent box
        for (var i = 0; i < image.numcomps; i++)
            _comps[i].bpcc = image.comps[i].prec - 1u + ((image.comps[i].sgnd ? 1u : 0u) << 7);

        //Color Specification box
        if (image.icc_profile_buf != null)
        {
            _meth = 2;
            _enumcs = 0;
        }
        else
        {
            _meth = 1;
            _enumcs = (uint)image.color_space;
        }

        uint alphaCount = 0, alphaChannel = 0, colorChannels = 0;
        for (uint i = 0; i < image.numcomps; i++)
        {
            if (image.comps[i].alpha != 0)
            {
                alphaCount++;
                alphaChannel = i;
            }
        }
        if (alphaCount == 1U)
        { /* no way to deal with more than 1 alpha channel */
            switch (_enumcs)
            {
                case 16:
                case 18:
                    colorChannels = 3;
                    break;
                case 17:
                    colorChannels = 1;
                    break;
                default:
                    alphaCount = 0U;
                    break;
            }
            if (alphaCount == 0U)
            {
                _cinfo.Warn("Alpha channel specified but unknown enumcs. No cdef box will be created.");
            }
            else if (image.numcomps < colorChannels + 1)
            {
                _cinfo.Warn("Alpha channel specified but not enough image components for an automatic cdef box creation.");
                alphaCount = 0U;
            }
            else if (alphaChannel < colorChannels)
            {
                _cinfo.Warn("Alpha channel position conflicts with color channel. No cdef box will be created.");
                alphaCount = 0U;
            }
        }
        else if (alphaCount > 1)
        {
            _cinfo.Warn("Multiple alpha channels specified. No cdef box will be created.");
        }
        if (alphaCount == 1U)
        { /* if here, we know what we can do */
            _color ??= new JP2Color();
            _color.channel_definitions = GC.AllocateUninitializedArray<JP2cdefInfo>((int)image.numcomps);

            uint i = 0;
            for (; i < colorChannels; i++)
            {
                _color.channel_definitions[i].cn = (ushort)i;
                _color.channel_definitions[i].typ = 0;
                _color.channel_definitions[i].asoc = (ushort)(i + 1U);
            }
            for (; i < image.numcomps; i++)
            {
                if (image.comps[i].alpha != 0)
                { /* we'll be here exactly once */
                    _color.channel_definitions[i].cn = (ushort) i;
                    _color.channel_definitions[i].typ = 1; // Opacity channel
                    _color.channel_definitions[i].asoc = 0; // Apply alpha channel to the whole image
                }
                else
                {
                    /* Unknown channel */
                    _color.channel_definitions[i].cn = (ushort) i;
                    _color.channel_definitions[i].typ = (ushort) 65535U;
                    _color.channel_definitions[i].asoc = (ushort) 65535U;
                }
            }
        }


        //C# PdfLib requires this.
        if (image.channel_definitions != null)
        {
            _color ??= new JP2Color();

            //Overwrites any channeld definition set above, this because
            //the definitions supplied by Pdflib are the correct definitions.
            _color.channel_definitions = image.channel_definitions;
        }

        _precedence = 0;
        _approx = 0;

        //JPIP is not supported by this impl.
        return true;
    }

    /// <summary>
    /// Configures for decoding
    /// </summary>
    /// <param name="cio">File data</param>
    /// <param name="parameters">Configuration</param>
    /// <remarks>
    /// 2.5
    /// </remarks>
    internal void SetupDecode(Cio cio, DecompressionParameters parameters)
    {
        _cio = cio;
        _j2K.SetupDecode(cio, parameters);

        _color = new JP2Color
        {
            ignore_pclr_cmap_cdef = (parameters.flags & DecompressionParameters.DPARAMETERS.IGNORE_PCLR_CMAP_CDEF_FLAG) != 0,
            // This is a C# addition.
            ignore_cmap = parameters.IgnoreColorLookupTable
        };
    }

    #endregion

    //2.5.1 - opj_jp2_read_header
    internal bool ReadHeader(out JPXImage image)
    {
        //Snip decoding validation (NOP)

        //Snip setup_header_reading, this creates a list over
        //what functions to call, which is then called in a for loop.
        //Here we call them directly.

        if (!ReadHeaderProcedure())
        {
            image = null;
            return false;
        }
        if (!HasJp2H)
        {
            _cinfo.Error("JP2H box missing. Required");
            image = null;
            return false;
        }
        if (!_hasIhdr)
        {
            _cinfo.Error("IHDR box missing. Required");
            image = null;
            return false;
        }

        var ret = _j2K.ReadHeader(out image);

        // Set Image Color Space
        if (image != null)
        {
            if (_enumcs is >= 11 and <= 18 && _enumcs != 15)
                image.color_space = (COLOR_SPACE)_enumcs;
            else if (_enumcs == 24)
                image.color_space = COLOR_SPACE.eYCC;
            else
                image.color_space = COLOR_SPACE.UNKNOWN;

            if (_color.icc_profile_buf != null)
            {
                image.icc_profile_buf = _color.icc_profile_buf;
                _color.icc_profile_buf = null;
            }
        }
        return ret;
    }

    //2.5 - opj_jp2_read_header_procedure
    private bool ReadHeaderProcedure()
    {
        while (ReadBoxhdr(out var box, out var nBytesRead))
        {
            //Codestream box
            if (box.type == JP2_Marker.JP2C)
            {
                if ((_state & JP2_STATE.HEADER) != 0)
                {
                    _state |= JP2_STATE.CODESTREAM;
                    return true;
                }
                else
                {
                    _cinfo.Error("Badly placed jpeg codestream\n");
                    return false;
                }
            }
            else if (box.length == 0)
            {
                _cinfo.Error("Cannot handle box of undefined sizes\n");
                return false;
            }
            else if (box.length < nBytesRead)
            {
                _cinfo.Error("invalid box size {0} ({1})\n", box.length, box.type.ToString());
                return false;
            }

            var handler = FindHandler(box.type);
            var currentDataSize = box.length - (uint)nBytesRead;
            if (currentDataSize > _cio.BytesLeft)
            {
                _cinfo.Error("Invalid box size {0} for box '{1}'. Need {2} bytes, {3} bytes remaining",
                    (int)box.length, box.type.ToString(), (int)currentDataSize, _cio.BytesLeft);
                return false;
            }

            if (handler != null)
            {
                if (!handler(box))
                    return false;
            }
            else
            {
                var hadlerMisplaced = ImgFindHandler(box.type);
                if (hadlerMisplaced != null)
                {
                    _cinfo.Warn("Found a misplaced {0} box outside jp2h box\n", box.type.ToString());
                    if ((_state & JP2_STATE.HEADER) != 0)
                    {
                        // Read anyway, we already have jp2h
                        box.data_length = currentDataSize;
                        if (!hadlerMisplaced(box))
                            return false;

                        continue;
                    }
                    else
                    {
                        _cinfo.Warn("JPEG2000 Header box not read yet, {0} box will be ignored\n", box.type.ToString());
                    }
                }
                    
                // Skip unkown boxes
                _cio.Skip(currentDataSize);
            }
        }

        return true;
    }

    /// <summary>
    /// Find handeler to use for interpeting a JP2 box
    /// </summary>
    /// <param name="type">Type of box</param>
    /// <returns>Handler or null if handler not found</returns>
    /// <remarks>
    /// 2.5
    /// The original implementation is more complex, but seeing as
    /// there are only 3 possible handlers, we keep this simple.
    /// </remarks>
    private Handeler FindHandler(JP2_Marker type)
    {
        switch (type)
        {
            case JP2_Marker.JP: return ReadJp;
            case JP2_Marker.FTYP: return ReadFtyp;
            case JP2_Marker.JP2H: return ReadJp2H;
        }
        return null;
    }

    /// <summary>
    /// Find handeler to use for interpeting a JP2 box
    /// </summary>
    /// <param name="type">Type of box</param>
    /// <returns>Handler or null if handler not found</returns>
    /// <remarks>
    /// 2.5
    /// The original implementation is more complex, but seeing as
    /// there are only 6 possible handlers, we keep this simple.
    /// </remarks>
    private Handeler ImgFindHandler(JP2_Marker type)
    {
        switch (type)
        {
            case JP2_Marker.IHDR: return ReadIhdr;
            case JP2_Marker.COLR: return ReadColr;
            case JP2_Marker.BPCC: return ReadBpcc;
            case JP2_Marker.PCLR: return ReadPclr;
            case JP2_Marker.CMAP: return ReadCmap;
            case JP2_Marker.CDEF: return ReadCdef;
        }
        return null;
    }

    //2.5 - opj_jp2_end_decompress
    internal bool EndDecompress()
    {
        if (_cio.BytesLeft > 8)
            ReadHeaderProcedure();

        return _j2K.EndDecompress();
    }

    private delegate bool Handeler(JP2Box box);

    //2.5.1 - opj_jp2_apply_color_postprocessing
    private bool ApplyColorPostprocessing(JPXImage image)
    {
        if (_j2K.NumcompsToDecode != 0)
        {
            // Bypass all JP2 component transforms
            return true;
        }

        if (!_color.ignore_pclr_cmap_cdef)
        {
            if (!CheckColor(image))
                return false;

            if (_color.jp2_pclr != null)
            {
                /* Part 1, I.5.3.4: Either both or none : */
                if (_color.jp2_pclr.cmap == null)
                    _color.jp2_pclr = null;
                else
                {
                    if (!_color.ignore_cmap)
                    {
                        if (!ApplyPclr(image, _color, _cinfo))
                            return false;
                    }
                    else
                    {
                        //Added for Pdflib
                        image.color_info = _color;
                    }
                }
            }

            // Apply the color space if needed
            if (_color.channel_definitions != null)
            {
                ApplyCdef(image, _color);
            }
        }

        return true;
    }

    //2.5.1 - opj_jp2_decode
    internal bool Decode(JPXImage image)
    {
        if (image == null) 
            return false;

        if (!_j2K.Decode(image))
        {
            _cinfo.Error("Failed to decode the codestream in the JP2 file");
            return false;
        }

        return ApplyColorPostprocessing(image);
    }

    /// <summary>
    /// Decodes a single tile in the image
    /// </summary>
    /// <remarks>2.5.1 - opj_jp2_get_tile</remarks>
    internal bool Decode(JPXImage image, uint tileNr)
    {
        if (image == null)
            return false;

        _cinfo.Warn("JP2 box which are after the codestream will not be read by this function.");

        if (!_j2K.Decode(image, tileNr))
        {
            _cinfo.Error("Failed to decode the codestream in the JP2 file");
            return false;
        }

        return ApplyColorPostprocessing(image);
    }

    //2.5 - opj_jp2_check_color
    private bool CheckColor(JPXImage image)
    {
        /* testcase 4149.pdf.SIGSEGV.cf7.3501 */
        if (_color.channel_definitions != null)
        {
            var info = _color.channel_definitions;
            var nrChannels = image.numcomps;
            ushort i;

            // cdef applies to cmap channels if any
            if (_color.jp2_pclr is { cmap: not null })
            {
                nrChannels = _color.jp2_pclr.nr_channels;
            }

            for (i = 0; i < info.Length; i++)
            {
                if (info[i].cn >= nrChannels)
                {
                    _cinfo.Error("Invalid component index {0} (>= {1})", info[i].cn, nrChannels);
                    return false;
                }
                if (info[i].asoc == 65535U)
                {
                    continue;
                }
                if (info[i].asoc > 0 && info[i].asoc - 1 >= nrChannels)
                {
                    _cinfo.Error("Invalid component index {0} (>= {1})", info[i].asoc - 1, nrChannels);
                    return false;
                }
            }

            // issue 397
            // ISO 15444-1 states that if cdef is present, it shall contain a complete list of channel definitions. */
            var n = (ushort)_color.channel_definitions.Length;
            while (nrChannels > 0)
            {
                for (i = 0; i < n; ++i)
                {
                    if (info[i].cn == nrChannels - 1U)
                    {
                        break;
                    }
                }
                if (i == n)
                {
                    _cinfo.Error("Incomplete channel definitions.");
                    return false;
                }
                --nrChannels;
            }
        }

        /* testcases 451.pdf.SIGSEGV.f4c.3723, 451.pdf.SIGSEGV.5b5.3723 and
           66ea31acbb0f23a2bbc91f64d69a03f5_signal_sigsegv_13937c0_7030_5725.pdf */
        if (_color.jp2_pclr is { cmap: not null })
        {
            var nrChannels = _color.jp2_pclr.nr_channels;
            var cmap = _color.jp2_pclr.cmap;
            var isSane = true;

            /* verify that all original components match an existing one */
            for (var i = 0; i < nrChannels; i++)
            {
                if (cmap[i].cmp >= image.numcomps)
                {
                    _cinfo.Error("Invalid component index {0} (>= {1}).", cmap[i].cmp, image.numcomps);
                    isSane = false;
                }
            }

            var pcolUsage = GC.AllocateUninitializedArray<bool>(nrChannels);
            if (pcolUsage == null)
            {
                _cinfo.Error("Unexpected OOM.");
                return false;
            }
                
            /* verify that no component is targeted more than once */
            for (var i = 0; i < nrChannels; i++)
            {
                var mtyp = cmap[i].mtyp;
                var pcol = cmap[i].pcol;
                if (mtyp != 0 && mtyp != 1)
                {
                    _cinfo.Error("Invalid value for cmap[{0}].mtyp = {1}.", i, mtyp);
                    isSane = false;
                } 
                else if (pcol >= nrChannels)
                {
                    _cinfo.Error("Invalid component/palette index for direct mapping {0}.", pcol);
                    isSane = false;
                }
                else if (pcolUsage[pcol] && mtyp == 1)
                {
                    _cinfo.Error("Component {0} is mapped twice.", pcol);
                    isSane = false;
                }
                else if (mtyp == 0 && pcol != 0)
                {
                    /* I.5.3.5 PCOL: If the value of the MTYP field for this channel is 0, then
                     * the value of this field shall be 0. */
                    _cinfo.Error("Direct use at #{0} however pcol={1}.", i, pcol);
                    isSane = false;
                }
                else if (mtyp == 1 && pcol != i)
                {
                    // OpenJPEG implementation limitation. See assert(i == pcol);
                    // in opj_jp2_apply_pclr() 
                    _cinfo.Error("Implementation limitation: for palette mapping, "+
                                 "pcol[{0}] should be equal to {1}, but is equal "+
                                 "to {2}.", i, i, pcol);
                    isSane = false;
                }
                else
                    pcolUsage[pcol] = true;
            }
            /* verify that all components are targeted at least once */
            for (var i = 0; i < nrChannels; i++)
            {
                if (!pcolUsage[i] && cmap[i].mtyp != 0)
                {
                    _cinfo.Error("Component {0} doesn't have a mapping.", i);
                    isSane = false;
                }
            }
            // Issue 235/447 weird cmap
            if (isSane && image.numcomps == 1U)
            {
                for (var i = 0; i < nrChannels; i++)
                {
                    if (!pcolUsage[i])
                    {
                        isSane = false;
                        _cinfo.Warn("Component mapping seems wrong. Trying to correct.");
                        break;
                    }
                }
                if (!isSane)
                {
                    isSane = true;
                    for (var i = 0; i < nrChannels; i++)
                    {
                        cmap[i].mtyp = 1;
                        cmap[i].pcol = (byte)i;
                    }
                }
            }

            if (!isSane)
                return false;
        }

        return true;
    }

    //2.5
    internal bool SetDecodeArea(JPXImage image, int startX, int startY, int endX, int endY)
    {
        return _j2K.SetDecodeArea(image, startX, startY, endX, endY);
    }

    /// <summary>
    /// Apply collected palette data
    /// </summary>
    /// <remarks>2.5 - opj_jp2_apply_pclr</remarks>
    internal static bool ApplyPclr(JPXImage image, JP2Color color, CompressionInfo cinfo)
    {
        uint max;
        ushort i, cmp, pcol;

        var channelSize = color.jp2_pclr.channel_size;
        var channelSign = color.jp2_pclr.channel_sign;
        var entries = color.jp2_pclr.entries;
        var cmap = color.jp2_pclr.cmap;
        var nrChannels = color.jp2_pclr.nr_channels;

        for(i = 0; i < nrChannels; ++i)
        {
            //Palette mapping
            cmp = cmap[i].cmp;
            if (image.comps[cmp].data == null)
            {
                if (cinfo != null)
                    cinfo.Error("image->comps[{0}].data == NULL in opj_jp2_apply_pclr()", i);
                return false;
            }
        }

        var oldComps = image.comps;
        var newComps = GC.AllocateUninitializedArray<ImageComp>(nrChannels);

        for(i = 0; i < nrChannels; ++i)
        {
            pcol = cmap[i].pcol; cmp = cmap[i].cmp;

            /* Direct use */
            if (cmap[i].mtyp == 0)
            {
                Debug.Assert(pcol == 0);
                newComps[i] = (ImageComp)oldComps[cmp].Clone();
            }
            else
            {
                Debug.Assert(pcol == i);
                newComps[pcol] = (ImageComp)oldComps[cmp].Clone();    
            }

            /* Palette mapping: */
            newComps[pcol].data = GC.AllocateUninitializedArray<int>((int)(oldComps[cmp].w * oldComps[cmp].h));
            newComps[pcol].prec = channelSize[i];
            newComps[pcol].sgnd = channelSign[i] != 0;
        }
        var topK = color.jp2_pclr.nr_entries - 1;

        for(i = 0; i < nrChannels; ++i)
        {
            /* Palette mapping: */
            cmp = cmap[i].cmp; 
            pcol = cmap[i].pcol;
            var src = oldComps[cmp].data; 
            max = newComps[i].w * newComps[i].h;

            /* Direct use: */
            int[] dst;
            uint j;
            if (cmap[i].mtyp == 0)
            {
                dst = newComps[i].data;
                for (j = 0; j < max; j++)
                    dst[j] = src[j];
            }
            else
            {
                dst = newComps[pcol].data;

                for (j = 0; j < max; ++j)
                {
                    /* The index */
                    int k;
                    if ((k = src[j]) < 0)
                        k = 0;
                    else if (k > topK)
                        k = topK;

                    /* The colour */
                    dst[j] = (int)entries[k * nrChannels + pcol];
                }
            }
        }
        max = image.numcomps;
        for (i = 0; i < max; i++)
        {
            if (oldComps[i].data != null)
                oldComps[i].data = null;
        }

        image.comps = newComps;
        image.numcomps = nrChannels;

        color.jp2_pclr = null;

        return true;
    }

    //2.5 - opj_jp2_apply_cdef
    private void ApplyCdef(JPXImage image, JP2Color color)
    {
        ushort typ;

        var info = color.channel_definitions;
        if (info == null) return;

        for(ushort i = 0; i < info.Length; ++i)
        {
            var asoc = info[i].asoc;
            var cn = info[i].cn;

            if (cn >= image.numcomps)
            {
                _cinfo.Warn("opj_jp2_apply_cdef: cn ={0}, numcomps ={1}",
                    cn, image.numcomps);
                continue;
            }
            if (asoc == 0 || asoc == 65535)
            {
                image.comps[cn].alpha = info[i].typ;
                continue;
            }

            var acn = (ushort) (asoc - 1);
            if (acn >= image.numcomps)
            {
                _cinfo.Warn("opj_jp2_apply_cdef: acn={0}, numcomps={1}", 
                    acn, image.numcomps);
                continue;
            }

            // Swap only if color channel
            if (cn != acn && info[i].typ == 0)
            {
                //C# Org impl does memcopies, but it is dealing with structs.
                (image.comps[cn], image.comps[acn]) = (image.comps[acn], image.comps[cn]);

                // Swap channels in following channel definitions, don't
                // bother with j <= i that are already processed
                for (var j = (ushort)(i + 1); j < info.Length; j++)
                {
                    if (info[j].cn == cn)
                        info[j].cn = acn;
                    else if (info[j].cn == acn)
                        info[j].cn = cn;
                    // asoc is related to color index. Do not update
                }
            }

            image.comps[cn].alpha = info[i].typ;
        }
	        
        color.channel_definitions = null;
    }

    /// <summary>
    /// Reads the Jpeg2000 file Header box - JP2 Header box (this box contains other boxes).
    /// </summary>
    /// <param name="box">The data contained in the file header box.</param>
    /// <returns>True if the JP2 Header box was successfully reconized</returns>
    /// <remarks>
    /// 2.5 - opj_jp2_read_jp2h
    /// </remarks>
    private bool ReadJp2H(JP2Box box) 
    {
        // Make sure the box is well placed
        if ((_state & JP2_STATE.FILE_TYPE) != JP2_STATE.FILE_TYPE)
        {
            _cinfo.Error("The  box must be the first box in the file.");
            return false;
        }

        _imgState = JP2_IMG_STATE.NONE;

        // iterate while there is data
        var headerSize = box.length - 8;
        while (headerSize > 0)
        {
            if (!ReadBoxhdr_char(out box, out var boxSize, (int)headerSize))
            {
                _cinfo.Error("Stream error while reading JP2 Header box");
                return false;
            }

            if (box.length > headerSize)
            {
                _cinfo.Error("Stream error while reading JP2 Header box: box length is inconsistent.");
                return false;
            }

            var handler = ImgFindHandler(box.type);
            box.data_length = box.length - (uint) boxSize;

            if (handler != null)
            {
                var pos = _cio.Pos;
                if (!handler(box))
                    return false;
                if (_cio.Pos - pos < box.data_length)
                {
                    //C# OpenJpeg 2.5 effectivly does this as it reads all the data
                    //   for a box before calling the handler. 
                    _cinfo.Warn("{0} box has {1} bytes of junk data",
                        box.type, box.data_length - (_cio.Pos - pos));
                    _cio.Skip((uint)(box.data_length - (_cio.Pos - pos)));
                }
            }
            else
            {
                _imgState |= JP2_IMG_STATE.UNKNOWN;
                _cio.Skip(box.data_length);
            }

            headerSize -= box.length;
        }

        if (!_hasIhdr)
        {
            _cinfo.Error("Stream error while reading JP2 Header box: no 'ihdr' box");
            return false;
        }

        _state |= JP2_STATE.HEADER;

        return true;
    }

    //2.5 - opj_jp2_read_cmap
    private bool ReadCmap(JP2Box box)
    {
        /* Need nr_channels: */
        if (_color.jp2_pclr == null)
        {
            _cinfo.Error("Need to read a PCLR box before the CMAP box.");
            return false;
        }

        /* Part 1, I.5.3.5: 'There shall be at most one Component Mapping box
         * inside a JP2 Header box' :
         */
        if (_color.jp2_pclr.cmap != null)
        {
            _cinfo.Error("Only one CMAP box is allowed.");
            return false;
        }

        var nrChannels = _color.jp2_pclr.nr_channels;
        if (box.data_length < nrChannels * 4)
        {
            _cinfo.Error("Insufficient data for CMAP box.");
            return false;
        }

        var cmap = new JP2cmap_comp[nrChannels];

        for(ushort i = 0; i < nrChannels; ++i)
        {
            cmap[i].cmp = _cio.ReadUShort();
            cmap[i].mtyp = _cio.ReadByte();
            cmap[i].pcol = _cio.ReadByte();
        }
        _color.jp2_pclr.cmap = cmap;

        return true;
    }

    //2.5 - opj_jp2_read_pclr
    private bool ReadPclr(JP2Box box)
    {
        var orgPos = _cio.Pos;

        /* Part 1, I.5.3.4: 'There shall be at most one Palette box inside
         * a JP2 Header box' :
         */
        if(_color.jp2_pclr != null) return false;

        if (box.data_length < 3)
            return false;

        var nrEntries = _cio.ReadUShort() /* NE */;
        if (nrEntries == 0 || nrEntries > 1024)
        {
            _cinfo.Error("Invalid PCLR box. Reports {0} entries", nrEntries);
            return false;
        }

        ushort nrChannels = _cio.ReadByte() /* NPC */;
        if (nrChannels == 0)
        {
            _cinfo.Error("Invalid PCLR box. Reports 0 palette columns");
            return false;
        }

        if (box.data_length < 3 + nrChannels)
            return false;

        var entries = new uint[nrChannels * nrEntries];
        var channelSize = new byte[nrChannels];
        var channelSign = new byte[nrChannels];

        var jp2Pclr = new JP2pclr
        {
            channel_sign = channelSign,
            channel_size = channelSize,
            entries = entries,
            nr_entries = nrEntries,
            nr_channels = nrChannels,
            cmap = null
        };

        _color.jp2_pclr = jp2Pclr;

        for(var i = 0; i < nrChannels; ++i)
        {
            var uc = _cio.ReadByte(); // Bi
            channelSize[i] = (byte) ((uc & 0x7f) + 1);
            channelSign[i] = (byte) ((uc & 0x80) == 0x80 ? 1 : 0);
        }

        for(int j = 0, k = 0; j < nrEntries; ++j)
        {
            for(var i = 0; i < nrChannels; ++i)
            {
                var bytesToRead = (uint)(channelSize[i] + 7) >> 3;

                //mem-b2ace68c-1381.jp2 triggers this condition. File decodes
                //fine without this check.
                if (box.data_length < _cio.Pos - orgPos + bytesToRead)
                    return false;

                /* Cji */
                entries[k++] = unchecked(_cio.Read(bytesToRead));
            }
        }

        return true;
    }

    /// <summary>
    /// Channel defenition box
    /// </summary>
    /// <remarks>
    /// 2.5 - opj_jp2_read_cdef
    /// This box defines what channels are alpha channels and such
    /// </remarks>
    private bool ReadCdef(JP2Box box)
    {
        /* Part 1, I.5.3.6: 'The shall be at most one Channel Definition box
         * inside a JP2 Header box.'
         */
        if(_color.channel_definitions != null) return false;

        if (box.data_length < 2)
        {
            _cinfo.Error("Insufficient data for CDEF box.");
            return false;
        }

        var n = _cio.ReadUShort();
        if (n == 0)
        {
            _cinfo.Error("Number of channel description is equal to zero in CDEF box.");
            return false;
        }

        if (box.data_length < 2 + n * 6)
        {
            _cinfo.Error("Insufficient data for CDEF box.");
            return false;
        }

        var info = new JP2cdefInfo[n];
        _color.channel_definitions = info;

        for(ushort i = 0; i < n; ++i)
        {
            info[i].cn = _cio.ReadUShort();
            info[i].typ = _cio.ReadUShort();
            info[i].asoc = _cio.ReadUShort();
        }

        return true;
    }

    //2.5 - opj_jp2_read_colr
    private bool ReadColr(JP2Box box) 
    {
        if (box.data_length < 3)
        {
            _cinfo.Error("Bad COLR header box (bad size)");
            return false;
        }

        /* Part 1, I.5.3.3 : 'A conforming JP2 reader shall ignore all Colour
         * Specification boxes after the first.'
         */
        if (_color.HasColor)
        {
            _cinfo.Info("A conforming JP2 reader shall ignore all Colour Specification boxes after the first, so we ignore this one.");
            _cio.Skip(box.data_length);
            return true;
        }

        _meth = _cio.ReadByte();
        _precedence = _cio.ReadByte();
        _approx = _cio.ReadByte();

        if (_meth == 1)
        {
            if (box.data_length < 7)
            {
                _cinfo.Error("Bad COLR header box (bad size: {0})", box.data_length);
                return false;
            }
            if (box.data_length > 7 && _enumcs != 14)
            {
                // Testcase Altona_Technical_v20_x4.pdf
                _cinfo.Warn("Bad COLR header box (bad size: {0})", box.data_length);
            }
            _enumcs = _cio.ReadUInt();

            if (_enumcs == 14)
            { // CIELab
                var cielab = new uint[9];
                cielab[0] = 14; // Enumcs

                uint ol, ra, oa, rb, ob;
                var rl = ra = rb = ol = oa = ob = 0;
                uint il = 0x00443530; // D50
                cielab[1] = 0x44454600; // DEF

                if (box.data_length == 35)
                {
                    rl = _cio.ReadUInt();
                    ol = _cio.ReadUInt();
                    ra = _cio.ReadUInt();
                    oa = _cio.ReadUInt();
                    rb = _cio.ReadUInt();
                    ob = _cio.ReadUInt();
                    il = _cio.ReadUInt();

                    cielab[1] = 0;
                } else if (box.data_length != 7)
                {
                    _cinfo.Warn("Bad COLR header box (CIELab, bad size: {0})", box.data_length);
                }
                cielab[2] = rl;
                cielab[4] = ra;
                cielab[6] = rb;
                cielab[3] = ol;
                cielab[5] = oa;
                cielab[7] = ob;
                cielab[8] = il;

                _color.icc_cielab_buf = cielab;
            }

            _color.HasColor = true;
        } 
        else if (_meth == 2)
        {
            /* ICC profile */
            var iccLen = (int) box.data_length - 3;
            Debug.Assert((int) (box.init_pos + box.length - _cio.Pos) == box.data_length - 3);

            _color.icc_profile_buf = new byte[iccLen];
            if (_cio.Read(_color.icc_profile_buf, 0, iccLen) != iccLen)
                throw new EndOfStreamException();

            _color.HasColor = true;
        }
        else
        {
            /*	ISO/IEC 15444-1:2004 (E), Table I.9 ­ Legal METH values:
                conforming JP2 reader shall ignore the entire Colour Specification box.*/
            _cio.Skip(box.data_length - 3);
        }

        return true;
    }

    //2.5 - opj_jp2_read_bpcc
    private bool ReadBpcc(JP2Box box)
    {
        if (_bpc != 255)
            _cinfo.Warn("A BPCC header box is available although BPC given by the IHDR box ({0}) indicate components bit depth is constant", _bpc);

        if (box.data_length != _numcomps)
        {
            _cinfo.Error("Bad BPCC header box (bad size)");
            return false;
        }

        for (var i = 0; i < _numcomps; i++)
        {
            _comps[i].bpcc = _cio.ReadByte();
        }

        return true;
    }

    //2.5 - opj_jp2_read_ihdr
    private bool ReadIhdr(JP2Box box)
    {
        if (_comps!= null)
        {
            _cinfo.Warn("Ignoring ihdr box. First ihdr box already read");
            return true;
        }

        if (box.data_length != 14)
        {
            _cinfo.Error("Bad image header box (bad size)");
            return false;
        }

        //Width and height
        _h = _cio.ReadUInt();
        _w = _cio.ReadUInt();

        _numcomps = _cio.ReadUShort();

        if (_h < 1 || _w < 1 || _numcomps < 1)
        {
            _cinfo.Error("Wrong values for: w{0}) h({1}) numcomps({2}) (ihdr)", _w, _h, _numcomps);
            return false;
        }
        if (_numcomps - 1U >= 16384U)
        {
            // Unsigned underflow is well defined: 1U <= jp2->numcomps <= 16384U
            _cinfo.Error("Invalid number of components (ihdr)");
            return false;

        }

        _comps = new JP2Comps[_numcomps];

        _bpc = _cio.ReadByte();

        _c = _cio.ReadByte();

        if (_c != 7)
        {
            _cinfo.Info("JP2 IHDR box: compression type indicate that the file is not a conforming JP2 file ({0}) ", _c);
        }

        _unkC = _cio.ReadBool();
        _ipr = _cio.ReadBool();

        _j2K.CP.AllowDifferentBitDepthSign = _bpc == 255;
        _j2K._ihdr_w = _w;
        _j2K._ihdr_h = _h;
        _hasIhdr = true;

        return true;
    }

    /// <summary>
    /// Reads a a FTYP box - File type box
    /// </summary>
    /// <param name="box">The data contained in the FTYP box</param>
    /// <returns>True if the FTYP box is valid</returns>
    /// <remarks>
    /// 2.5 - opj_jp2_read_ftyp
    /// </remarks>
    private bool ReadFtyp(JP2Box box)
    {
        if (_state != JP2_STATE.SIGNATURE)
        {
            _cinfo.Error("The ftyp box must be the second box in the file.");
            return false;
        }

        if (box.length < 16)
        {
            _cinfo.Error("Error with FTYP signature Box size");
            return false;
        }

        _brand = (JP2_Marker)_cio.ReadUInt();
        _minversion = _cio.ReadUInt();

        var remainingBytes = (int) box.length - 16;

        // Number of bytes must be a multiple of 4
        if ((remainingBytes & 0x3) != 0)
        {
            _cinfo.Error("Error with FTYP signature Box size");
            return false;
        }

        _cl = new JP2_Marker[remainingBytes / 4];

        for (var i = 0; i < _cl.Length; i++)
        {
            _cl[i] = (JP2_Marker)_cio.ReadUInt();
        }

        _state |= JP2_STATE.FILE_TYPE;

        return true;
    }

    /// <summary>
    /// Reads a jpeg2000 file signature box.
    /// </summary>
    /// <param name="box">The data contained in the signature box</param>
    /// <returns>Rrue if the file signature box is valid</returns>
    /// <remarks>
    /// 2.5 - opj_jp2_read_jp
    /// </remarks>
    private bool ReadJp(JP2Box box)
    {
        if (_state != JP2_STATE.NONE)
        {
            _cinfo.Error("The signature box must be the first box in the file.");
            return false;
        }

        if (box.length != 12)
        {
            _cinfo.Error("Error with JP signature Box size");
            return false;
        }

        if (0x0d0a870a != _cio.ReadInt())
        {
            _cinfo.Error("Error with JP Signature : bad magic number");
            return false;
        }

        _state |= JP2_STATE.SIGNATURE;

        return true;
    }

    /// <summary>
    /// Reads a box header. The box is the way data is packed inside a jpeg2000 file structure.
    /// </summary>
    /// <remarks>2.5 - opj_jp2_read_boxhdr</remarks>
    private bool ReadBoxhdr(out JP2Box box, out int nBytesRead)
    {
        box = new JP2Box
        {
            init_pos = _cio.Pos
        };

        if (_cio.BytesLeft < 8)
        {
            nBytesRead = (int)_cio.BytesLeft;
            _cio.Skip((uint)nBytesRead);
            return false;
        }

        box.length = _cio.ReadUInt();
        box.type = (JP2_Marker) _cio.ReadUInt();
        nBytesRead = 8;

        // Do we have a "special very large box ?
        // read then the XLBo
        if (box.length == 1)
        {
            if (_cio.ReadInt() != 0)
            {
                _cinfo.Error("Cannot handle box sizes higher than 2^32");
                nBytesRead += 4;
                return false;
            }
            box.length = _cio.ReadUInt();
            nBytesRead = 16;
        }
        else if (box.length == 0) // last box
        {
            var bleft = _cio.BytesLeft;
            if (bleft > 0xFFFFFFFFL - 8L)
            {
                _cinfo.Error("Cannot handle box sizes higher than 2^32");
                return false;
            }
            box.length = (uint) (bleft + 8);
        }

        return true;
    }

    //2.5 - opj_jp2_read_boxhdr_char
    private bool ReadBoxhdr_char(out JP2Box box, out int nBytesRead, int maxSize)
    {
        box = new JP2Box();
        if (maxSize < 8)
        {
            nBytesRead = 0;
            _cinfo.Error("Cannot handle box of less than 8 bytes");
            return false;
        }

        box.init_pos = _cio.Pos;
        box.length = _cio.ReadUInt();
        box.type = (JP2_Marker)_cio.ReadUInt();
        nBytesRead = 8;

        // Do we have a "special very large box
        // read then the XLBox
        if (box.length == 1)
        {
            if (maxSize < 8)
            {
                _cinfo.Error("Cannot handle XL box of less than 16 bytes");
                return false;
            }

            if (_cio.ReadInt() != 0)
            {
                _cinfo.Error("Cannot handle box sizes higher than 2^32");
                nBytesRead += 4;
                return false;
            }
            box.length = _cio.ReadUInt();
            nBytesRead = 16;

            if (box.length == 0)
            {
                _cinfo.Error("Cannot handle box of undefined sizes");
                return false;
            }
        }
        else if (box.length == 0)
        {
            _cinfo.Error("Cannot handle box of undefined sizes");
            return false;
        }
        if (box.length < nBytesRead)
        {
            _cinfo.Error("Box length is inconsistent");
            return false;
        }

        return true;
    }

    //2.5 - opj_jp2_encode
    internal bool Encode()
    {
        return _j2K.Encode();
    }

    //2.5 - opj_jp2_default_validation
    private bool DefaultValidation(Stream cio)
    {
        var lIsValid = true;

        /* JPEG2000 codec validation */

        /* STATE checking */
        /* make sure the state is at 0 */
        lIsValid &= _state == JP2_STATE.NONE;

        /* make sure not reading a jp2h ???? WEIRD */
        lIsValid &= _imgState == JP2_IMG_STATE.NONE;

        /* POINTER validation */
        /* make sure a j2k codec is present */
        lIsValid &= _j2K != null;

        /* make sure a procedure list is present */
        //l_is_valid &= (_procedure_list != null);

        /* make sure a validation list is present */
        //l_is_valid &= (_validation_list != null);

        /* PARAMETER VALIDATION */
        /* number of components */
        lIsValid &= _cl != null;
        /* width */
        lIsValid &= _h > 0;
        /* height */
        lIsValid &= _w > 0;
        /* precision */
        for (var i = 0; i < _numcomps; ++i)
        {
            lIsValid &= (_comps[i].bpcc & 0x7FU) < 38U; //Bug in org. impl?
        }

        /* METH */
        lIsValid &= _meth is > 0 and < 3;

        /* stream validation */
        /* back and forth is needed */
        lIsValid &= cio.CanSeek;

        return lIsValid;
    }

    //2.5 - opj_jp2_start_compress
    internal bool StartCompress(Cio cio)
    {
        var bcio = new BufferCio(cio);
        {
            var buf = new byte[256];
            bcio.SetBuffer(ref buf, 256);
        }

        if (!DefaultValidation(cio.Stream))
            return false;

        if (!WriteHeader(bcio))
            return false;

        //Makes room for the Code Stream marker, which will be
        //written later.
        Debug.Assert(bcio.BufferPos == 0);
        SkipJp2C(cio.Stream);

        return _j2K.StartCompress(bcio);
    }

    //2.5 - opj_jp2_end_compress
    internal bool EndCompress()
    {
        var bcio = _j2K.EndGetBCIO();

        // Writes header
        WriteJp2C(bcio);

        return true;
    }

    //2.5 - opj_jp2_setup_header_writing
    private bool WriteHeader(BufferCio bcio)
    {
        WriteJp(bcio);
        WriteFtyp(bcio);
        if (!WriteJp2H(bcio))
            return false;

        //C# Skip is called by the parent function
        return true;
    }

    /// <summary>
    /// Makes room for the Code Stream marker and
    /// store away the position.
    /// </summary>
    /// <remarks>2.5 - opj_jp2_skip_jp2c</remarks>
    private void SkipJp2C(Stream cio)
    {
        _j2KCodestreamOffset = cio.Position;
        cio.Seek(8, SeekOrigin.Current);
    }


    /// <summary>
    /// Writes the Jpeg2000 codestream Header box - JP2C Header box. This function must be called AFTER the coding has been done.
    /// </summary>
    /// <remarks>2.5 - opj_jp2_write_jp2c</remarks>
    private void WriteJp2C(BufferCio bcio) 
    {
        /* J2K encoding */
        var j2KCodestreamExit = bcio.Pos;
        var j2KCodestreamLength = j2KCodestreamExit - _j2KCodestreamOffset;

        bcio.Pos = _j2KCodestreamOffset;
        bcio.Write((int) j2KCodestreamLength);
        bcio.Write(JP2_Marker.JP2C);
        bcio.Commit();
        bcio.Pos = j2KCodestreamExit;
    }

    /// <summary>
    /// JP2 header
    /// </summary>
    /// <remarks>2.5 - opj_jp2_write_jp2h</remarks>
    private bool WriteJp2H(BufferCio bcio)
    {
        uint jp2HSize = 8;

        //C# - Calculate the needed buffer size
        //IHDR needs 22 bytes
        jp2HSize += 22;

        //Adds for Bit per Component
        uint bpccSize = 0;
        if (_bpc == 255)
            bpccSize += 8 + _numcomps;
        jp2HSize += bpccSize;

        //Adds for the color info
        uint colrSize = 11;
        switch(_meth)
        {
            case 1:
                colrSize += 4; break;
            case 2:
                colrSize += (uint) _color.icc_profile_buf.Length; break;
            default:
                return false;
        }
        jp2HSize += colrSize;

        uint cdefSize = 0;
        if (_color is { channel_definitions: not null })
            cdefSize = 10u + 6u * (uint)_color.channel_definitions.Length;
        jp2HSize += cdefSize;

        bcio.SetBuffer(jp2HSize);

        //Writes out the length
        bcio.Write(jp2HSize);

        //Signature
        bcio.Write(JP2_Marker.JP2H);

        WriteIhdr(bcio);

        if (_bpc == 255)
            Write_BPCC(bcio, bpccSize);
        WriteColr(bcio, colrSize);

        if (_color is { channel_definitions: not null })
            WriteCdef(bcio, cdefSize);

        Debug.Assert(bcio.BufferPos == jp2HSize);
        bcio.Commit();

        return true;
    }

    /// <summary>
    ///  Writes the Channel Definition box.
    /// </summary>
    /// <remarks>2.5 - opj_jp2_write_cdef</remarks>
    private void WriteCdef(BufferCio bcio, uint cdefSize)
    {
        var channelDefinitions = _color.channel_definitions;

        bcio.Write(cdefSize);
        bcio.Write(JP2_Marker.CDEF);

        //Writes number of definitions
        bcio.WriteUShort(channelDefinitions.Length);

        for (var c = 0; c < channelDefinitions.Length; c++)
        {
            bcio.WriteUShort(channelDefinitions[c].cn);
            bcio.WriteUShort(channelDefinitions[c].typ);
            bcio.WriteUShort(channelDefinitions[c].asoc);
        }
    }

    /// <summary>
    /// Writes the Image Header box - Image Header box.
    /// </summary>
    /// <remarks>2.5 - opj_jp2_write_ihdr</remarks>
    private void WriteIhdr(BufferCio bcio)
    {
        bcio.Write(22); // Size of the box
        bcio.Write(JP2_Marker.IHDR);

        //Writes out the height and width
        bcio.Write(_h);
        bcio.Write(_w);

        bcio.Write(_numcomps, 2);

        bcio.Write(_bpc, 1);	

        bcio.Write(_c, 1);
        bcio.Write(_unkC);
        bcio.Write(_ipr);
    }

    /// <summary>
    /// Writes the Bit per Component box
    /// </summary>
    /// <remarks>2.5 - opj_jp2_write_bpcc</remarks>
    private void Write_BPCC(BufferCio bcio, uint bpccSize)
    {
        bcio.Write(bpccSize);
        bcio.Write(JP2_Marker.BPCC);

        for (var i = 0; i < _numcomps; i++)
            bcio.Write(_comps[i].bpcc, 1);
    }

    /// <summary>
    /// Writes the Colour Specification box
    /// </summary>
    /// <remarks>2.5 - opj_jp2_write_colr</remarks>
    private void WriteColr(BufferCio bcio, uint colrSize)
    {
        bcio.Write(colrSize);
        bcio.Write(JP2_Marker.COLR);

        bcio.Write(_meth, 1);
        bcio.Write(_precedence, 1);
        bcio.Write(_approx, 1);

        if (_meth == 1)
            bcio.Write(_enumcs);
        else
            bcio.Write(_color.icc_profile_buf, 0, _color.icc_profile_buf.Length);
    }

    /// <summary>
    /// File type
    /// </summary>
    /// <remarks>2.5 - opj_jp2_write_ftyp</remarks>
    private void WriteFtyp(BufferCio bcio)
    {
        var ftypSize = 16 + 4 * _cl.Length;
        bcio.SetBuffer((uint)ftypSize);

        //Writes the length
        bcio.Write(ftypSize);

        //Signature
        bcio.Write(JP2_Marker.FTYP);

        bcio.Write(_brand);
        bcio.Write(_minversion);

        for (var i = 0; i < _cl.Length; i++)
            bcio.Write(_cl[i]);

        bcio.Commit();
    }

    /// <summary>
    /// Signature
    /// </summary>
    /// <remarks>2.5 - opj_jp2_write_jp</remarks>
    private bool WriteJp(BufferCio bcio)
    {
        //Writes the length
        bcio.Write(12);

        //Writes out the JP2 signature
        bcio.Write(JP2_Marker.JP);

        //Writes magic number
        bcio.Write(0x0d0a870a);

        bcio.Commit();

        return true;
    }
}