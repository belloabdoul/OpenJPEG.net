﻿#region License

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

namespace OpenJpeg;

/// <summary>
/// This class contains all information needed
/// for compression or decompression
/// </summary>
public sealed class CompressionInfo
{
    /// <summary>
    /// Version number of this libary
    /// </summary>
    public static string Version => "2.5.3";

    /// <summary>
    /// Whenever this object is a decompressor or
    /// a compressor
    /// </summary>
    private readonly bool _isDecompressor;

    /// <summary>
    /// Type of compression or decompression
    /// </summary>
    private CodecFormat _codecFormat;

    /// <summary>
    /// Handles j2k stuff
    /// </summary>
    private readonly J2K _j2K;

    /// <summary>
    /// Handles JP2 stuff
    /// </summary>
    private readonly Jp2 _jp2;

    /// <summary>
    /// For sending messages back to the client
    /// </summary>
    private EventMgr _mgr;

    /// <summary>
    /// Get or set event manager
    /// </summary>
    public EventMgr EventManager
    {
        get => _mgr;
        set => _mgr = value;
    }

    /// <summary>
    /// Optional data object for events
    /// </summary>
    public object ClientData { get; }

    /// <summary>
    /// Whenever this object is set up for decompression
    /// </summary>
    public bool IsDecompressor => _isDecompressor;

    /// <summary>
    /// Allow multithreading
    /// </summary>
    internal bool DisableMultiThreading { get; set; }

    /// <summary>
    /// Functions for compressing an image
    /// </summary>
    /// <remarks>Openjpeg 2.1 API</remarks>
    private readonly Compression _cFuncs;

    /// <summary>
    /// Functions for decompresssing an image
    /// </summary>
    /// <remarks>Openjpeg 2.1 API</remarks>
    private readonly Decompression _dFuncs;

    /// <summary>
    /// Constructor
    /// </summary>
    /// <param name="isDecompressor">
    /// Set true if you intent to decompress
    /// </param>
    /// <remarks>
    /// opj_create_compress
    /// </remarks>
    public CompressionInfo(bool isDecompressor, CodecFormat format)
    {
        _isDecompressor = isDecompressor;
        _codecFormat = format;
        _j2K = J2K.Create(this);
        if (format == CodecFormat.Jpeg2P)
        {
            _jp2 = new Jp2(this, _j2K);
            if (_isDecompressor)
            {
                _dFuncs.SetupDecoder = _jp2.SetupDecode;
                _dFuncs.ReadHeader = _jp2.ReadHeader;
                _dFuncs.SetDecodeArea = _jp2.SetDecodeArea;
                _dFuncs.Decode = _jp2.Decode;
                _dFuncs.TileDecode = _jp2.Decode;
                _dFuncs.EndDecompress = _jp2.EndDecompress;
            }
            else
            {
                _cFuncs.SetupEncoder = _jp2.SetupEncoder;
                _cFuncs.StartCompress = _jp2.StartCompress;
                _cFuncs.EndCompress = _jp2.EndCompress;
                _cFuncs.Encode = _jp2.Encode;
            }
        }
        else if (format == CodecFormat.Jpeg2K)
        {
            //_j2k = J2K.Create(this);
            if (_isDecompressor)
            {
                _dFuncs.SetupDecoder = _j2K.SetupDecode;
                _dFuncs.ReadHeader = _j2K.ReadHeader;
                _dFuncs.SetDecodeArea = _j2K.SetDecodeArea;
                _dFuncs.Decode = _j2K.Decode;
                _dFuncs.TileDecode = _j2K.Decode;
                _dFuncs.EndDecompress = _j2K.EndDecompress;
            }
            else
            {
                _cFuncs.SetupEncoder = _j2K.SetupEncoder;
                _cFuncs.StartCompress = _j2K.StartCompress;
                _cFuncs.EndCompress = _j2K.EndCompress;
                _cFuncs.Encode = _j2K.Encode;
            }
        }
        else
            throw new NotSupportedException(format.ToString());
    }

    /// <summary>
    /// Configures the decoder
    /// </summary>
    /// <param name="cio">File data</param>
    /// <param name="parameters">Configuration</param>
    /// <returns>True if the setup was a success</returns>
    /// <remarks>
    /// 2.5
    /// </remarks>
    public bool SetupDecoder(Cio cio, DecompressionParameters parameters)
    {
        if (cio is { CanSeek: true } && parameters != null && _isDecompressor)
        {
            DisableMultiThreading = parameters.DisableMultiThreading;
            _dFuncs.SetupDecoder(cio, parameters);
            return true;
        }

        return false;
    }

    /// <summary>
    /// Reads the header of the file
    /// </summary>
    /// <param name="image">Image with header information</param>
    /// <returns>True on success</returns>
    /// <remarks>
    /// 2.5
    /// </remarks>
    public bool ReadHeader(out JPXImage image)
    {
        if (_isDecompressor)
        {
            return _dFuncs.ReadHeader(out image);
        }
        else
        {
            image = null;
            return false;
        }
    }

    //2.5
    public bool SetDecodeArea(JPXImage image, int startX, int startY, int endX, int endY)
    {
        if (_isDecompressor)
            return _dFuncs.SetDecodeArea(image, startX, startY, endX, endY);
        return false;
    }

    /// <summary>
    /// Sets up the encoder. Required before encoding.
    /// </summary>
    /// <param name="parameters"></param>
    /// <param name="image"></param>
    /// <returns>True if setup was a sucess</returns>
    public bool SetupEncoder(CompressionParameters parameters, JPXImage image)
    {
        if (parameters != null && image != null && !_isDecompressor && parameters.Valid)
        {
            DisableMultiThreading = parameters.DisableMultiThreading;
            return _cFuncs.SetupEncoder(parameters, image);
        }

        return false;
    }

    public bool SetExtraOptions(ExtraOption extra)
    {
        if (extra != null && _j2K != null)
            return _j2K.SetExtraOptions(extra);

        return false;
    }

    //2.5 - opj_start_compress
    public bool StartCompress(Cio cio)
    {
        if (cio != null)
        {
            if (!_isDecompressor)
                return _cFuncs.StartCompress(cio);
        }

        return false;
    }

    public bool Encode()
    {
        if (!_isDecompressor)
            return _cFuncs.Encode();
        return false;
    }

    /// <summary>
    /// End to compress the current image
    /// </summary>
    /// <remarks>2.5 - opj_end_compress</remarks>
    public bool EndCompress()
    {
        if (!_isDecompressor)
            return _cFuncs.EndCompress();
        return false;
    }

    /// <summary>
    /// Opens a CIO stream
    /// </summary>
    /// <param name="file">File to read or write</param>
    /// <param name="read">Whenever to read or write</param>
    /// <remarks>
    /// Openjpeg 2.1 now supports streams. Yay. But since my impl.
    /// already supports streams, I'll keep this old 1.4 method.
    /// 
    /// I.e. this is now eq. with "opj_stream_create_default_file_stream"
    /// </remarks>
    public Cio OpenCio(Stream file, bool read)
    {
        ArgumentNullException.ThrowIfNull(file);

        return read ? new Cio(this, file, OpenMode.Read) : new Cio(this, file, _j2K.ImageLength);
    }

    //2.5
    public bool Decode(JPXImage image)
    {
        if (image == null) throw new ArgumentNullException();

        if (_isDecompressor)
            return _dFuncs.Decode(image);

        return false;
    }

    public bool Decode(JPXImage image, uint tileNo)
    {
        if (image == null) throw new ArgumentNullException();

        if (_isDecompressor)
            return _dFuncs.TileDecode(image, tileNo);

        return false;
    }

    //2.5 - opj_end_decompress
    public bool EndDecompress()
    {
        if (_isDecompressor)
        {
            return _dFuncs.EndDecompress();
        }

        return false;
    }

    /// <summary>
    /// Emits a formated error string
    /// </summary>
    /// <param name="msg">Message</param>
    /// <param name="arg0">Parameters to put into the message</param>
    internal void Error(string msg, params object[] arg0)
    {
        if (_mgr == null || _mgr.Error == null) return;

        _mgr.Error(string.Format(msg, arg0), ClientData);
    }

    /// <summary>
    /// Emits a formated error string
    /// </summary>
    /// <param name="msg">Message</param>
    /// <param name="arg0">Parameters to put into the message</param>
    internal void ErrorMt(string msg, params object[] arg0)
    {
        if (_mgr == null || _mgr.Error == null) return;

        //Locking on _mgr is not ideal, as this object is publically visible
        lock (_mgr) _mgr.Error(string.Format(msg, arg0), ClientData);
    }

    /// <summary>
    /// Emits a formated info string
    /// </summary>
    /// <param name="msg">Message</param>
    /// <param name="arg0">Parameters to put into the message</param>
    internal void Info(string msg, params object[] arg0)
    {
        if (_mgr == null || _mgr.Info == null) return;

        _mgr.Info(string.Format(msg, arg0), ClientData);
    }

    /// <summary>
    /// Emits a formated warning string
    /// </summary>
    /// <param name="msg">Message</param>
    /// <param name="arg0">Parameters to put into the message</param>
    internal void Warn(string msg, params object[] arg0)
    {
        if (_mgr == null || _mgr.Warning == null) return;

        _mgr.Warning(string.Format(msg, arg0), ClientData);
    }

    /// <summary>
    /// Emits a formated warning string
    /// </summary>
    /// <param name="msg">Message</param>
    /// <param name="arg0">Parameters to put into the message</param>
    internal void WarnMt(string msg, params object[] arg0)
    {
        if (_mgr == null || _mgr.Warning == null) return;

        //Locking on _mgr is not ideal, as this object is publically visible
        lock (_mgr) _mgr.Warning(string.Format(msg, arg0), ClientData);
    }
}

/// <summary>
/// OpenJpeg 2.1 collects the relevant compression functions into
/// a struct.
/// </summary>
internal struct Compression
{
    internal delegate bool StartCompressFunc(Cio cio);

    internal delegate bool EncodeFunc();

    internal delegate bool WriteTileFunc(uint tileIndex, byte[] data, int dataSize, Cio cio);

    internal delegate bool EndCompressFunc();

    internal delegate void DestroyFunc();

    internal delegate bool SetupEncoderFunc(CompressionParameters cparameters, JPXImage image);

    internal StartCompressFunc StartCompress;
    internal EncodeFunc Encode;
    internal WriteTileFunc WriteTile;
    internal EndCompressFunc EndCompress;
    internal DestroyFunc Destroy;
    internal SetupEncoderFunc SetupEncoder;
}

internal struct Decompression
{
    internal delegate void SetupDecoderFunc(Cio cio, DecompressionParameters param);

    internal delegate bool ReadHeaderFunc(out JPXImage image);

    internal delegate bool SetDecodeAreaFunc(JPXImage image, int startX, int startY, int endX, int endY);

    internal delegate bool DecodeFunc(JPXImage image);

    internal delegate bool TileDecodeFunc(JPXImage image, uint tileIndex);

    internal delegate bool EndDecompressFunc();

    internal SetupDecoderFunc SetupDecoder;
    internal ReadHeaderFunc ReadHeader;
    internal SetDecodeAreaFunc SetDecodeArea;
    internal DecodeFunc Decode;
    internal TileDecodeFunc TileDecode;
    internal EndDecompressFunc EndDecompress;
}