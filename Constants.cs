namespace OpenJpeg;

public static class Constants
{
    /// <summary>
    /// Number of maximum resolution level authorized
    /// </summary>
    internal const int J2KMaxrlvls = 33;

    /// <summary>
    /// Number of maximum sub-band linked to number of resolution level
    /// </summary>
    internal const int J2KMaxbands = 3 * J2KMaxrlvls - 2;
    internal const int J2KDefaultNbSegs = 10;
    internal const int J2KTcdMatrixMaxLayerCount = 10;
    internal const int J2KTcdMatrixMaxResolutionCount = 10;

    internal const int MqcNumctxs = 19;

    internal const int T1Numctxs = (int)T1_CTXNO.UNI + (int)T1_NUMCTXS.UNI;

    internal const int T1NmsedecBits = 7;

    internal const int T1NmsedecFracbits = T1NmsedecBits - 1;

    public const int Cinema24Cs = 1302083;
    public const int Cinema48Cs = 651041;
    public const int Cinema24Comp = 1041666;
    public const int Cinema48Comp = 520833;

    internal const int MctDefaultNbRecords = 10;
    internal const int MccDefaultNbRecords = 10;

    internal const int DefaultCblkDataSize = 8192;

    /// <summary>
    /// Org impl have this as ulong.MaxValue, but C# has a 2GB limit.
    /// </summary>
    internal const int SizeMax = int.MaxValue;

    /// <summary>
    /// Default size of the buffer for header bytes
    /// </summary>
    internal const int J2KDefaultHeaderSize = 1000;

    /// <summary>
    /// Margin for a fake FFFF marker
    /// </summary>
    internal const int CommonCblkDataExtra = 2;

    internal const int ParamDefaultNumresolution = 6;
    internal const int ParamDefaultCblockw = 64;
    internal const int ParamDefaultCblockh = 64;
    internal const PROG_ORDER ParamDefaultProgOrder = PROG_ORDER.LRCP;

    /// <summary>
    /// Maximum main level
    /// </summary>
    internal const ushort ImfMainlevelMax = 11;
}