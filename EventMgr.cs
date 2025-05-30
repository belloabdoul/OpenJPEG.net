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
using System.Diagnostics;

namespace OpenJpeg;

/// <summary>
/// Implement this delegate to consume msg events
/// </summary>
/// <param name="msg">The message</param>
/// <param name="clientData">Optional client data</param>
public delegate void MsgCallback(string msg, object clientData);

public class EventMgr
{
    #region Variables and properties

    internal readonly MsgCallback Error;
    internal readonly MsgCallback Warning;
    internal readonly MsgCallback Info;

    #endregion

    #region Init

    /// <summary>
    /// Constructor
    /// </summary>
    /// <param name="error">Error message callback</param>
    /// <param name="warn">Warning message callback</param>
    /// <param name="info">Information message callback</param>
    public EventMgr(MsgCallback error, MsgCallback warn, MsgCallback info)
    { Error = error; Warning = warn; Info = info; }

    #endregion

    /// <summary>
    /// Open Jpeg clock.
    /// </summary>
    /// <remarks>
    /// Put here for convinience. 
    /// </remarks>
    internal static double obj_clock()
    {
        return Stopwatch.GetTimestamp() / (double)Stopwatch.Frequency;
    }
}