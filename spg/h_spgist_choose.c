/*-----Choose function*/
Datum
spg_kmer_choose(PG_FUNCTION_ARGS)
{
    spgChooseIn *in = (spgChooseIn *) PG_GETARG_POINTER(0);
    spgChooseOut *out = (spgChooseOut *) PG_GETARG_POINTER(1);

    kmer *inKmer = DatumGetKmerP(in->datum);
    char *inData = inKmer->data;
    int inK = inKmer->k;
    char *prefixData = NULL;
    int prefixK = 0;
    int commonLen = 0;
    int16 nodeChar = 0;
    int i = 0;

    /* Check for prefix match, set nodeChar to first byte after prefix */
    if (in->hasPrefix)
    {
        kmer *prefixKmer = DatumGetKmerP(in->prefixDatum);
        prefixData = prefixKmer->data;
        prefixK = prefixKmer->k;
        
        commonLen = commonPrefix(inData + in->level, prefixData, inK - in->level, prefixK);

        if (commonLen == prefixK)
        {
            /* node label --- first non-common character */
            if (inK - in->level > commonLen)
                nodeChar = *(unsigned char *) (inData + in->level + commonLen);
            else
                nodeChar = -1; /*completely common values */
        }
        else /*Not match prefix  -> split tuple*/
        {            
            out->resultType = spgSplitTuple;

            if (commonLen == 0)
            {
                out->result.splitTuple.prefixHasPrefix = false;
            }
            else
            {
                out->result.splitTuple.prefixHasPrefix = true;
                out->result.splitTuple.prefixPrefixDatum =
                    formKmerDatum(prefixData, commonLen);
            }
            out->result.splitTuple.prefixNNodes = 1;
            out->result.splitTuple.prefixNodeLabels =
                (Datum *) palloc(sizeof(Datum));
            out->result.splitTuple.prefixNodeLabels[0] =
                Int16GetDatum(*(unsigned char *) (prefixData + commonLen));

            out->result.splitTuple.childNodeN = 0;

            if (prefixK - commonLen == 1)
            {
                out->result.splitTuple.postfixHasPrefix = false;
            }
            else
            {
                out->result.splitTuple.postfixHasPrefix = true;
                out->result.splitTuple.postfixPrefixDatum =
                    formKmerDatum(prefixData + commonLen + 1, prefixK - commonLen - 1);
            }

            PG_RETURN_VOID();
        }
    }
    else if (inK > in->level)
    {        
        nodeChar = *(unsigned char *) (inData + in->level); /* node label = 1st character after the current level */
    }
    else
    {
        nodeChar = -1;
    }

    /* Look up nodeChar in the node label array */
    if (searchChar(in->nodeLabels, in->nNodes, nodeChar, &i))
    {
        /* Descend to the existing node */
        int levelAdd;

        out->resultType = spgMatchNode;
        out->result.matchNode.nodeN = i;
        levelAdd = commonLen;
        if (nodeChar >= 0)
            levelAdd++;
        out->result.matchNode.levelAdd = levelAdd;
        if (inK - in->level - levelAdd > 0)
        {
            out->result.matchNode.restDatum =
                formKmerDatum(inData + in->level + levelAdd,
                              inK - in->level - levelAdd);
        }
        else
        {
            out->result.matchNode.restDatum = formKmerDatum(NULL, 0);
        }
    }
    else if (in->allTheSame)
    {
        /* Cannot use AddNode; split the tuple */
        out->resultType = spgSplitTuple;
        out->result.splitTuple.prefixHasPrefix = in->hasPrefix;
        out->result.splitTuple.prefixPrefixDatum = in->prefixDatum;
        out->result.splitTuple.prefixNNodes = 1;
        out->result.splitTuple.prefixNodeLabels = (Datum *) palloc(sizeof(Datum));
        out->result.splitTuple.prefixNodeLabels[0] = Int16GetDatum(-2);
        out->result.splitTuple.childNodeN = 0;
        out->result.splitTuple.postfixHasPrefix = false;
    }
    else
    {
        /* Add a new node for the unseen nodeChar */
        out->resultType = spgAddNode;
        out->result.addNode.nodeLabel = Int16GetDatum(nodeChar);
        out->result.addNode.nodeN = i;
    }

    PG_RETURN_VOID();
}
