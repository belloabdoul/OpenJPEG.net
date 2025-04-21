using System;

namespace OpenJpeg.Internal;

internal sealed class SparseArrayInt32
{
    public readonly uint Width;
    public readonly uint Height;
    public readonly uint BlockWidth;
    public readonly uint BlockHeight;
    public readonly uint BlockCountHor;
    public readonly uint BlockCountVer;
    public readonly int[][] DataBlocks;

    private SparseArrayInt32(uint w, uint h, uint bw, uint bh, uint bch, uint bcv,
        int[][] db)
    {
        Width = w;
        Height = h;
        BlockWidth = bw;
        BlockHeight = bh;
        BlockCountHor = bch;
        BlockCountVer = bcv;
        DataBlocks= db;
    }

    private bool is_region_valid(uint x0,
        uint y0,
        uint x1,
        uint y1)
    {
        return !(x0 >= Width || x1 <= x0 || x1 > Width ||
                 y0 >= Height || y1 <= y0 || y1 > Height);
    }

    /// <remarks>
    /// 2.5 - opj_sparse_array_int32_read_or_write
    /// 
    /// C# Uses a float for a buffer, but treats it like an int[]
    /// </remarks>
    private bool read_or_write(uint x0,
        uint y0,
        uint x1,
        uint y1,
        float[] buf,
        int bufPos,
        uint bufColStride,
        uint bufLineStride,
        bool forgiving,
        bool isReadOp)
    {
        uint yIncr = 0;
        var iof = new IntOrFloat();

        if (!is_region_valid(x0, y0, x1, y1))
        {
            return forgiving;
        }

        var blockY = y0 / BlockHeight;
        for (var y = y0; y < y1; blockY++, y += yIncr)
        {
            uint xIncr = 0;
            yIncr = y == y0 ? BlockHeight - y0 % BlockHeight :
                BlockHeight;
            var blockYOffset = BlockHeight - yIncr;
            yIncr = Math.Min(yIncr, y1 - y);
            var blockX = x0 / BlockWidth;
            for (var x = x0; x < x1; blockX++, x += xIncr)
            {
                uint j;
                xIncr = x == x0 ? BlockWidth - x0 % BlockWidth : BlockWidth;
                var blockXOffset = BlockWidth - xIncr;
                xIncr = Math.Min(xIncr, x1 - x);
                var srcBlock = DataBlocks[blockY * BlockCountHor + blockX];
                if (isReadOp)
                {
                    if (srcBlock == null)
                    {
                        if (bufColStride == 1)
                        {
                            var destPtr = (int)(bufPos + (y - y0) * bufLineStride +
                                                 (x - x0) * bufColStride);
                            for (j = 0; j < yIncr; j++)
                            {
                                Array.Clear(buf, destPtr, (int)xIncr);
                                destPtr += (int)bufLineStride;
                            }
                        }
                        else
                        {
                            var destPtr = (int)(bufPos + (y - y0) * bufLineStride +
                                                 (x - x0) * bufColStride);
                            for (j = 0; j < yIncr; j++)
                            {
                                for (uint k = 0; k < xIncr; k++)
                                {
                                    buf[destPtr + k * bufColStride] = 0;
                                }
                                destPtr += (int)bufLineStride;
                            }
                        }
                    }
                    else
                    {
                        var srcPtr = blockYOffset * BlockWidth + blockXOffset;
                        if (bufColStride == 1)
                        {
                            var destPtr = (int)(bufPos + (y - y0) * bufLineStride
                                                         +
                                                         (x - x0) * bufColStride);
                            //C# Commented out since it dosn't make sense for .net
                            //if (x_incr == 4)
                            //{
                            //    /* Same code as general branch, but the compiler */
                            //    /* can have an efficient memcpy() */
                            //    (void)(x_incr); /* trick to silent cppcheck duplicateBranch warning */
                            //    for (j = 0; j < y_incr; j++)
                            //    {
                            //        memcpy(dest_ptr, src_ptr, sizeof(OPJ_INT32) * x_incr);
                            //        dest_ptr += buf_line_stride;
                            //        src_ptr += block_width;
                            //    }
                            //}
                            //else
                            {
                                for (j = 0; j < yIncr; j++)
                                {
                                    Buffer.BlockCopy(srcBlock, (int)srcPtr * sizeof(int), buf, destPtr * sizeof(int), sizeof(int) * (int)xIncr);
                                    destPtr += (int)bufLineStride;
                                    srcPtr += BlockWidth;
                                }
                            }
                        }
                        else
                        {
                            var destPtr = (uint)bufPos + (y - y0) * bufLineStride
                                                         +
                                                         (x - x0) * bufColStride;
                            if (xIncr == 1)
                            {
                                for (j = 0; j < yIncr; j++)
                                {
                                    buf[destPtr] = srcBlock[srcPtr];
                                    destPtr += bufLineStride;
                                    srcPtr += BlockWidth;
                                }
                            }
                            else if (yIncr == 1 && bufColStride == 2)
                            {
                                uint k;
                                for (k = 0; k < (xIncr & ~3U); k += 4)
                                {
                                    buf[destPtr + k * bufColStride] = srcBlock[srcPtr + k];
                                    buf[destPtr + (k + 1) * bufColStride] = srcBlock[srcPtr + k + 1];
                                    buf[destPtr + (k + 2) * bufColStride] = srcBlock[srcPtr + k + 2];
                                    buf[destPtr + (k + 3) * bufColStride] = srcBlock[srcPtr + k + 3];
                                }
                                for (; k < xIncr; k++)
                                {
                                    buf[destPtr + k * bufColStride] = srcBlock[srcPtr + k];
                                }
                            }
                            else if (xIncr >= 8 && bufColStride == 8)
                            {
                                for (j = 0; j < yIncr; j++)
                                {
                                    uint k;
                                    for (k = 0; k < (xIncr & ~3U); k += 4)
                                    {
                                        buf[destPtr + k * bufColStride] = srcBlock[srcPtr + k];
                                        buf[destPtr + (k + 1) * bufColStride] = srcBlock[srcPtr + k + 1];
                                        buf[destPtr + (k + 2) * bufColStride] = srcBlock[srcPtr + k + 2];
                                        buf[destPtr + (k + 3) * bufColStride] = srcBlock[srcPtr + k + 3];
                                    }
                                    for (; k < xIncr; k++)
                                    {
                                        buf[destPtr + k * bufColStride] = srcBlock[srcPtr + k];
                                    }
                                    destPtr += bufLineStride;
                                    srcPtr += BlockWidth;
                                }
                            }
                            else
                            {
                                /* General case */
                                for (j = 0; j < yIncr; j++)
                                {
                                    for (uint k = 0; k < xIncr; k++)
                                    {
                                        buf[destPtr + k * bufColStride] = srcBlock[srcPtr + k];
                                    }
                                    destPtr += bufLineStride;
                                    srcPtr += BlockWidth;
                                }
                            }
                        }
                    }
                }
                else
                {
                    if (srcBlock == null)
                    {
                        srcBlock = new int[BlockWidth * BlockHeight];
                        DataBlocks[blockY * BlockCountHor + blockX] = srcBlock;
                    }

                    if (bufColStride == 1)
                    {
                        var destPtr = (int)(0 + blockYOffset *
                            BlockWidth + blockXOffset);
                        var srcPtr = (int)(bufPos + (y - y0) *
                            bufLineStride + (x - x0) * bufColStride);
                        //if (x_incr == 4)
                        //{
                        //    /* Same code as general branch, but the compiler */
                        //    /* can have an efficient memcpy() */
                        //    (void)(x_incr); /* trick to silent cppcheck duplicateBranch warning */
                        //    for (j = 0; j < y_incr; j++)
                        //    {
                        //        memcpy(dest_ptr, src_ptr, sizeof(OPJ_INT32) * x_incr);
                        //        dest_ptr += block_width;
                        //        src_ptr += buf_line_stride;
                        //    }
                        //}
                        //else
                        {
                            for (j = 0; j < yIncr; j++)
                            {
                                Buffer.BlockCopy(buf, srcPtr * sizeof(int), srcBlock, destPtr * sizeof(int), (int)xIncr * sizeof(int));
                                destPtr += (int)BlockWidth;
                                srcPtr += (int)bufLineStride;
                            }
                        }
                    }
                    else
                    {
                        var destPtr = 0 + blockYOffset *
                            BlockWidth + blockXOffset;
                        var srcPtr = (uint)bufPos + (y - y0) *
                            bufLineStride + (x - x0) * bufColStride;
                        if (xIncr == 1)
                        {
                            for (j = 0; j < yIncr; j++)
                            {
                                iof.F = buf[srcPtr];
                                srcBlock[destPtr] = iof.I;
                                srcPtr += bufLineStride;
                                destPtr += BlockWidth;
                            }
                        }
                        else if (xIncr >= 8 && bufColStride == 8)
                        {
                            for (j = 0; j < yIncr; j++)
                            {
                                uint k;
                                for (k = 0; k < (xIncr & ~3U); k += 4)
                                {
                                    iof.F = buf[srcPtr + k * bufColStride];
                                    srcBlock[destPtr + k] = iof.I;
                                    iof.F = buf[srcPtr + (k + 1) * bufColStride];
                                    srcBlock[destPtr + k + 1] = iof.I;
                                    iof.F = buf[srcPtr + (k + 2) * bufColStride];
                                    srcBlock[destPtr + k + 2] = iof.I;
                                    iof.F = buf[srcPtr + (k + 3) * bufColStride];
                                    srcBlock[destPtr + k + 3] = iof.I;
                                }
                                for (; k < xIncr; k++)
                                {
                                    iof.F = buf[srcPtr + k * bufColStride];
                                    srcBlock[destPtr + k] = iof.I;
                                }
                                srcPtr += bufLineStride;
                                destPtr += BlockWidth;
                            }
                        }
                        else
                        {
                            /* General case */
                            for (j = 0; j < yIncr; j++)
                            {
                                for (uint k = 0; k < xIncr; k++)
                                {
                                    iof.F = buf[srcPtr + k * bufColStride];
                                    srcBlock[destPtr + k] = iof.I;
                                }
                                srcPtr += bufLineStride;
                                destPtr += BlockWidth;
                            }
                        }
                    }
                }
            }
        }

        return true;
    }

    //2.5 - opj_sparse_array_int32_read_or_write
    private bool read_or_write(uint x0,
        uint y0,
        uint x1,
        uint y1,
        int[] buf,
        int bufPos,
        uint bufColStride,
        uint bufLineStride,
        bool forgiving,
        bool isReadOp)
    {
        uint yIncr;

        if (!is_region_valid(x0, y0, x1, y1))
        {
            return forgiving;
        }

        var blockY = y0 / BlockHeight;
        for (var y = y0; y < y1; blockY++, y += yIncr)
        {
            uint xIncr;
            yIncr = y == y0 ? BlockHeight - y0 % BlockHeight :
                BlockHeight;
            var blockYOffset = BlockHeight - yIncr;
            yIncr = Math.Min(yIncr, y1 - y);
            var blockX = x0 / BlockWidth;
            for (var x = x0; x < x1; blockX++, x += xIncr)
            {
                uint j;
                xIncr = x == x0 ? BlockWidth - x0 % BlockWidth : BlockWidth;
                var blockXOffset = BlockWidth - xIncr;
                xIncr = Math.Min(xIncr, x1 - x);
                var srcBlock = DataBlocks[blockY * BlockCountHor + blockX];
                if (isReadOp)
                {
                    if (srcBlock == null)
                    {
                        if (bufColStride == 1)
                        {
                            var destPtr = (int)(bufPos + (y - y0) * bufLineStride +
                                                 (x - x0) * bufColStride);
                            for (j = 0; j < yIncr; j++)
                            {
                                Array.Clear(buf, destPtr, (int)xIncr);
                                destPtr += (int)bufLineStride;
                            }
                        }
                        else
                        {
                            var destPtr = (int)(bufPos + (y - y0) * bufLineStride +
                                                 (x - x0) * bufColStride);
                            for (j = 0; j < yIncr; j++)
                            {
                                for (uint k = 0; k < xIncr; k++)
                                {
                                    buf[destPtr + k * bufColStride] = 0;
                                }
                                destPtr += (int)bufLineStride;
                            }
                        }
                    }
                    else
                    {
                        var srcPtr = blockYOffset * BlockWidth + blockXOffset;
                        if (bufColStride == 1)
                        {
                            var destPtr = (int)(bufPos + (y - y0) * bufLineStride
                                                         +
                                                         (x - x0) * bufColStride);
                            //if (x_incr == 4)
                            //{
                            //    /* Same code as general branch, but the compiler */
                            //    /* can have an efficient memcpy() */
                            //    (void)(x_incr); /* trick to silent cppcheck duplicateBranch warning */
                            //    for (j = 0; j < y_incr; j++)
                            //    {
                            //        memcpy(dest_ptr, src_ptr, sizeof(OPJ_INT32) * x_incr);
                            //        dest_ptr += buf_line_stride;
                            //        src_ptr += block_width;
                            //    }
                            //}
                            //else
                            {
                                for (j = 0; j < yIncr; j++)
                                {
                                    Buffer.BlockCopy(srcBlock, (int)srcPtr * sizeof(int), buf, destPtr * sizeof(int), sizeof(int) * (int)xIncr);
                                    destPtr += (int)bufLineStride;
                                    srcPtr += BlockWidth;
                                }
                            }
                        }
                        else
                        {
                            var destPtr = (uint)bufPos + (y - y0) * bufLineStride
                                                         +
                                                         (x - x0) * bufColStride;
                            if (xIncr == 1)
                            {
                                for (j = 0; j < yIncr; j++)
                                {
                                    buf[destPtr] = srcBlock[srcPtr];
                                    destPtr += bufLineStride;
                                    srcPtr += BlockWidth;
                                }
                            }
                            else if (yIncr == 1 && bufColStride == 2)
                            {
                                uint k;
                                for (k = 0; k < (xIncr & ~3U); k += 4)
                                {
                                    buf[destPtr + k * bufColStride] = srcBlock[srcPtr + k];
                                    buf[destPtr + (k + 1) * bufColStride] = srcBlock[srcPtr + k + 1];
                                    buf[destPtr + (k + 2) * bufColStride] = srcBlock[srcPtr + k + 2];
                                    buf[destPtr + (k + 3) * bufColStride] = srcBlock[srcPtr + k + 3];
                                }
                                for (; k < xIncr; k++)
                                {
                                    buf[destPtr + k * bufColStride] = srcBlock[srcPtr + k];
                                }
                            }
                            else if (xIncr >= 8 && bufColStride == 8)
                            {
                                for (j = 0; j < yIncr; j++)
                                {
                                    uint k;
                                    for (k = 0; k < (xIncr & ~3U); k += 4)
                                    {
                                        buf[destPtr + k * bufColStride] = srcBlock[srcPtr + k];
                                        buf[destPtr + (k + 1) * bufColStride] = srcBlock[srcPtr + k + 1];
                                        buf[destPtr + (k + 2) * bufColStride] = srcBlock[srcPtr + k + 2];
                                        buf[destPtr + (k + 3) * bufColStride] = srcBlock[srcPtr + k + 3];
                                    }
                                    for (; k < xIncr; k++)
                                    {
                                        buf[destPtr + k * bufColStride] = srcBlock[srcPtr + k];
                                    }
                                    destPtr += bufLineStride;
                                    srcPtr += BlockWidth;
                                }
                            }
                            else
                            {
                                /* General case */
                                for (j = 0; j < yIncr; j++)
                                {
                                    for (uint k = 0; k < xIncr; k++)
                                    {
                                        buf[destPtr + k * bufColStride] = srcBlock[srcPtr + k];
                                    }
                                    destPtr += bufLineStride;
                                    srcPtr += BlockWidth;
                                }
                            }
                        }
                    }
                }
                else
                {
                    if (srcBlock == null)
                    {
                        srcBlock = GC.AllocateUninitializedArray<int>((int)(BlockWidth * BlockHeight));
                        DataBlocks[blockY * BlockCountHor + blockX] = srcBlock;
                    }

                    if (bufColStride == 1)
                    {
                        var destPtr = (int)(0 + blockYOffset *
                            BlockWidth + blockXOffset);
                        var srcPtr = (int)(bufPos + (y - y0) *
                            bufLineStride + (x - x0) * bufColStride);
                        //if (x_incr == 4)
                        //{
                        //    /* Same code as general branch, but the compiler */
                        //    /* can have an efficient memcpy() */
                        //    (void)(x_incr); /* trick to silent cppcheck duplicateBranch warning */
                        //    for (j = 0; j < y_incr; j++)
                        //    {
                        //        memcpy(dest_ptr, src_ptr, sizeof(OPJ_INT32) * x_incr);
                        //        dest_ptr += block_width;
                        //        src_ptr += buf_line_stride;
                        //    }
                        //}
                        //else
                        {
                            for (j = 0; j < yIncr; j++)
                            {
                                Buffer.BlockCopy(buf, srcPtr * sizeof(int), srcBlock, destPtr * sizeof(int), (int)xIncr * sizeof(int));
                                destPtr += (int)BlockWidth;
                                srcPtr += (int)bufLineStride;
                            }
                        }
                    }
                    else
                    {
                        var destPtr = 0 + blockYOffset *
                            BlockWidth + blockXOffset;
                        var srcPtr = (uint)bufPos + (y - y0) *
                            bufLineStride + (x - x0) * bufColStride;
                        if (xIncr == 1)
                        {
                            for (j = 0; j < yIncr; j++)
                            {
                                srcBlock[destPtr] = buf[srcPtr];
                                srcPtr += bufLineStride;
                                destPtr += BlockWidth;
                            }
                        }
                        else if (xIncr >= 8 && bufColStride == 8)
                        {
                            for (j = 0; j < yIncr; j++)
                            {
                                uint k;
                                for (k = 0; k < (xIncr & ~3U); k += 4)
                                {
                                    srcBlock[destPtr + k] = buf[srcPtr + k * bufColStride];
                                    srcBlock[destPtr + k + 1] = buf[srcPtr + (k + 1) * bufColStride];
                                    srcBlock[destPtr + k + 2] = buf[srcPtr + (k + 2) * bufColStride];
                                    srcBlock[destPtr + k + 3] = buf[srcPtr + (k + 3) * bufColStride];
                                }
                                for (; k < xIncr; k++)
                                {
                                    srcBlock[destPtr + k] = buf[srcPtr + k * bufColStride];
                                }
                                srcPtr += bufLineStride;
                                destPtr += BlockWidth;
                            }
                        }
                        else
                        {
                            /* General case */
                            for (j = 0; j < yIncr; j++)
                            {
                                for (uint k = 0; k < xIncr; k++)
                                {
                                    srcBlock[destPtr + k] = buf[srcPtr + k * bufColStride];
                                }
                                srcPtr += bufLineStride;
                                destPtr += BlockWidth;
                            }
                        }
                    }
                }
            }
        }

        return true;
    }
    internal bool Read(uint x0,
        uint y0,
        uint x1,
        uint y1,
        int[] dest,
        int destPos,
        uint destColStride,
        uint destLineStride,
        bool forgiving)
    {
        return read_or_write(x0, y0, x1, y1,
            dest,
            destPos,
            destColStride,
            destLineStride,
            forgiving,
            true);
    }

    //2.5 - opj_sparse_array_int32_read
    internal bool Read(uint x0,
        uint y0,
        uint x1,
        uint y1,
        float[] dest,
        int destPos,
        uint destColStride,
        uint destLineStride,
        bool forgiving)
    {
        return read_or_write(x0, y0, x1, y1,
            dest,
            destPos,
            destColStride,
            destLineStride,
            forgiving,
            true);
    }

    internal bool Write(uint x0,
        uint y0,
        uint x1,
        uint y1,
        int[] src,
        int srcPos,
        uint srcColStride,
        uint srcLineStride,
        bool forgiving)
    {
        return read_or_write(x0, y0, x1, y1,
            src,
            srcPos,
            srcColStride,
            srcLineStride,
            forgiving,
            false);
    }

    internal bool Write(uint x0,
        uint y0,
        uint x1,
        uint y1,
        float[] src,
        int srcPos,
        uint srcColStride,
        uint srcLineStride,
        bool forgiving)
    {
        return read_or_write(x0, y0, x1, y1,
            src,
            srcPos,
            srcColStride,
            srcLineStride,
            forgiving,
            false);
    }

    internal static SparseArrayInt32 Create(uint width, uint height, 
        uint blockWidth, uint blockHeight)
    {
        if (width == 0 || height == 0 || blockWidth == 0 || blockHeight == 0)
        {
            return null;
        }
        if (blockWidth > ~0U / blockHeight / sizeof(int))
        {
            return null;
        }

        var bch = MyMath.uint_ceildiv(width, blockWidth);
        var bcv = MyMath.uint_ceildiv(height, blockHeight);
        if (bch > ~0U / bcv)
        {
            return null;
        }

        return new SparseArrayInt32(width, height, blockWidth, blockHeight,
            bch, bcv, new int[bch * bcv][]);
    }

    internal static SparseArrayInt32 Init(TcdTilecomp tilec, uint numres)
    {
        var trMax = tilec.resolutions[numres - 1];
        var w = (uint)(trMax.x1 - trMax.x0);
        var h = (uint)(trMax.y1 - trMax.y0);

        var sa = Create(w, h, Math.Min(w, 64u), Math.Min(h, 64u));
        if (sa == null)
            return null;

        for (uint resno = 0; resno < numres; ++resno)
        {
            var res = tilec.resolutions[resno];

            for (uint bandno = 0; bandno < res.numbands; ++bandno)
            {
                var band = res.bands[bandno];

                for (uint precno = 0; precno < res.pw * res.ph; ++precno)
                {
                    var precinct = band.precincts[precno];
                    for (uint cblkno = 0; cblkno < precinct.cw * precinct.ch; ++cblkno)
                    {
                        var cblk = precinct.dec[cblkno];
                        if (cblk.decoded_data != null)
                        {
                            var x = (uint)(cblk.x0 - band.x0);
                            var y = (uint)(cblk.y0 - band.y0);
                            var cblkW = (uint)(cblk.x1 - cblk.x0);
                            var cblkH = (uint)(cblk.y1 - cblk.y0);

                            if ((band.bandno & 1) != 0)
                            {
                                var pres = tilec.resolutions[resno - 1];
                                x += (uint)(pres.x1 - pres.x0);
                            }
                            if ((band.bandno & 2) != 0)
                            {
                                var pres = tilec.resolutions[resno - 1];
                                y += (uint)(pres.y1 - pres.y0);
                            }

                            if (!sa.Write(x, y,
                                    x + cblkW, y + cblkH,
                                    cblk.decoded_data, 0,
                                    1, cblkW, true))
                            {
                                return null;
                            }
                        }
                    }
                }
            }
        }

        return sa;
    }
}