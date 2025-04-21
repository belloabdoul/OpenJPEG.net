﻿using System;
using System.IO;

namespace OpenJpeg.Util;

//Based on WBIO just without the "skip bit" feature
internal class BitWriter
{
    internal byte[] _buf = GC.AllocateUninitializedArray<byte>(256);
    private int _buf_pos;
    private byte _unfinished_byte;
    private int _u_pos = 8;
    private readonly Stream _target;

    internal BitWriter(Stream target)
    {
        _target = target;
    }

    public void Write(int value, int nbits)
    {
        for (var i = nbits - 1; i >= 0; i--)
            WriteBit((value >> i) & 1);
    }

    public void WriteBit(int bit)
    {
        //Commits the buffer when full
        if (_u_pos == 0)
        {
            _u_pos = 8;

            //Writes out the bytes
            _buf[_buf_pos++] = _unfinished_byte;
            if (_buf_pos == _buf.Length)
            {
                _target.Write(_buf, 0, _buf_pos);
                _buf_pos = 0;
            }
            _unfinished_byte = 0;
        }
            
        //Appends a bit
        _unfinished_byte |= (byte) (bit << --_u_pos);
    }
    public void WriteBit(bool bit) { WriteBit(bit ? 1 : 0); }

    /**
     * Empties the buffers
     *
     * Data is always byte aligned after flushing.
     *
     * @return False if not all data could be written out.
     */
    public bool Flush()
    {
        if (_u_pos < 8)
            _buf[_buf_pos++] = _unfinished_byte;
        _u_pos = 8; _unfinished_byte = 0;
        if (_buf_pos > 0)
        {
            _target.Write(_buf, 0, _buf_pos);
            _buf_pos = 0;
        }

        return true;
    }
}