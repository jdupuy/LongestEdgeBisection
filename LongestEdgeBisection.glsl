/* leb.glsl - public domain Longest Edge Bisection GLSL library
by Jonathan Dupuy

*/

#ifndef BUFFER_BINDING_LEB
#   error User must specify the binding of the LEB buffer
#endif
layout(std430, binding = BUFFER_BINDING_LEB)
buffer LebBuffer {
    uint u_LebBuffer[];
};

#ifndef LEB_MAX_DEPTH
#   error User must specify the maximum LoD for the subdivision
#endif

#define LEB_BUFFER u_LebBuffer

// data structures
struct leb_Node {
    uint id;    // binary code
    int depth;  // subdivision depth
};
struct leb_SameDepthNeighborIDs {
    uint left, right, edge, _reserved;
};
struct leb_DiamondParent {
    leb_Node base, top;
};
struct leb_NodeAndNeighbors {
    leb_Node left, right, edge, node;
};

// manipulation
void leb_SplitNodeConforming(in const leb_Node node);
void leb_MergeNodeConforming(in const leb_Node node,
                             in const leb_DiamondParent diamond);

// O(1) queries
uint leb_NodeCount();
bool leb_IsLeafNode(in const leb_Node node);
bool leb_IsRootNode(in const leb_Node node);
bool leb_IsNullNode(in const leb_Node node);
leb_Node leb_ParentNode(in const leb_Node node);
leb_SameDepthNeighborIDs leb_GetSameDepthNeighborIDs(in const leb_NodeAndNeighbors nodes);

// O(depth) queries
uint                     leb_EncodeNode(in const leb_Node node);
leb_Node                 leb_DecodeNode(uint nodeID);
leb_NodeAndNeighbors     leb_DecodeNodeAndNeighbors(uint nodeID);
leb_SameDepthNeighborIDs leb_DecodeSameDepthNeighborIDs(in const leb_Node node);
leb_DiamondParent        leb_DecodeDiamondParent(in const leb_Node node);

// subdivision routine O(depth)
vec3   leb_DecodeNodeAttributeArray(in const leb_Node node, in const vec3 data);
mat2x3 leb_DecodeNodeAttributeArray(in const leb_Node node, in const mat2x3 data);
mat3x3 leb_DecodeNodeAttributeArray(in const leb_Node node, in const mat3x3 data);
mat4x3 leb_DecodeNodeAttributeArray(in const leb_Node node, in const mat4x3 data);

// intersection test O(depth)
leb_Node leb_BoundingNode(vec2 p);

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------


/*******************************************************************************
 * GetBitValue -- Returns the value of a bit stored in a 32-bit word
 *
 */
uint leb__GetBitValue(uint bitField, uint bitID)
{
    return ((bitField >> bitID) & 1u);
}


/*******************************************************************************
 * SetBitValue -- Sets the value of a bit stored in a 32-bit word
 *
 */
void leb__SetBitValue(uint bufferIndex, uint bitID, uint bitValue)
{
    const uint bitMask = ~(1u << bitID);

    atomicAnd(LEB_BUFFER[bufferIndex], bitMask);
    atomicOr(LEB_BUFFER[bufferIndex], bitValue << bitID);
}


/*******************************************************************************
 * BitFieldInsert -- Returns the bit field after insertion of some bit data in range
 * [bitOffset, bitOffset + bitCount - 1]
 *
 */
void leb__BitFieldInsert(uint bufferIndex, uint bitData, uint bitOffset, uint bitCount)
{
    uint bitMask = ~(~(0xFFFFFFFFu << bitCount) << bitOffset);

    atomicAnd(LEB_BUFFER[bufferIndex], bitMask);
    atomicOr(LEB_BUFFER[bufferIndex], bitData << bitOffset);
}


/*******************************************************************************
 * BitFieldExtract -- Extracts bits [bitOffset, bitOffset + bitCount - 1] from
 * a bit field, returning them in the least significant bits of the result.
 *
 */
uint leb__BitFieldExtract(uint bitField, uint bitOffset, uint bitCount)
{
    uint bitMask = ~(0xFFFFFFFFu << bitCount);

    return (bitField >> bitOffset) & bitMask;
}


/*******************************************************************************
 * NodeToLeafBitID -- Computes the bitfield index associated to a node
 *
 */
uint leb__NodeToLeafBitID(in const leb_Node node)
{
    return (node.id << (LEB_MAX_DEPTH - node.depth)) + (2u << LEB_MAX_DEPTH);
}


/*******************************************************************************
 * SetNodeBufferValue -- Sets the bit associated to a leaf node to bitValue
 *
 */
void leb__SetNodeBitValue(in const leb_Node node, uint bitValue)
{
    uint bitID = leb__NodeToLeafBitID(node);
    uint bufferID = bitID >> 5u;
    uint localBitID = bitID & 31u;

    leb__SetBitValue(bufferID, localBitID, bitValue);
}


/*******************************************************************************
 * GetNodeBufferValue -- Returns the value of the bit associated to a leaf node
 *
 */
uint leb__GetNodeBitValue(in const leb_Node node)
{
    uint bitID = leb__NodeToLeafBitID(node);
    uint bufferID = bitID >> 5u;
    uint bitOffset = bitID & 31u;

    return leb__GetBitValue(LEB_BUFFER[bufferID], bitOffset);
}


/*******************************************************************************
 * BitFieldSize -- Computes the number of uints to allocate for the buffer
 *
 * For a tree of max depth D, the number of bits is 2^(D+2).
 * Note that 2 bits are "wasted" in the sense that they only serve
 * to round the required number of bits to a power of two.
 *
 */
uint leb__BitFieldSize(uint lebMaxDepth)
{
    return 1u << (lebMaxDepth + 2u);
}


/*******************************************************************************
 * BufferUintSize -- Computes the number of uints to allocate for the bitfield
 *
 */
uint leb__BufferUintSize(uint lebMaxDepth)
{
    return leb__BitFieldSize(lebMaxDepth) >> 5u;
}


/*******************************************************************************
 * DataBitID -- Returns the bitID associated to dataID
 *
 * For a LEB of max depth D and given an index in [0, 2^(D+1) - 1], this
 * functions is used to emulate the behaviour of a lookup in an array, i.e.,
 * uint[dataID]. It provides the first bit in the bitfield that stores
 * information associated with the element of index dataID.
 *
 * For data located at level d, the bit offset is 2^d x (3 - d + D)
 * We then offset this quantity by the index by (dataID - 2^d) x (D + 1 - d)
 * Note that the null index (dataID = 0) is also supported.
 *
 */
uint leb__DataBitID(uint dataID, int dataDepth)
{
    int tmp = 1 << dataDepth;
    int bitCount = 1 + LEB_MAX_DEPTH - dataDepth;
    int offset = tmp * (2 + bitCount);
    int elementID = int(dataID) - tmp;

    return uint(offset + elementID * bitCount);
}


/*******************************************************************************
 * DataBitSize -- Returns the number of bits associated with a given node
 *
 */
int leb__DataBitSize(int nodeDepth)
{
    return LEB_MAX_DEPTH - nodeDepth + 1;
}


/*******************************************************************************
 * SetData -- Sets bitCount bits located at arrayID to bitData
 *
 * Note that this procedure writes to at most two uint32 elements.
 * Two elements are relevant whenever the specified interval overflows 32-bit
 * words.
 *
 */
void
leb__SetDataExplicit(uint nodeID, int nodeDepth, uint bitCount, uint bitData)
{
    uint alignedBitOffset = leb__DataBitID(nodeID, nodeDepth);
    uint maxBufferIndex = leb__BufferUintSize(LEB_MAX_DEPTH) - 1u;
    uint bufferIndexLSB = (alignedBitOffset >> 5u);
    uint bufferIndexMSB = min(bufferIndexLSB + 1, maxBufferIndex);
    uint bitFieldOffsetLSB = alignedBitOffset & 31u;
    uint bitCountLSB = min(32u - bitFieldOffsetLSB, bitCount);
    uint bitCountMSB = bitCount - bitCountLSB;

    leb__BitFieldInsert(bufferIndexLSB, bitData               , bitFieldOffsetLSB, bitCountLSB);
    leb__BitFieldInsert(bufferIndexMSB, bitData >> bitCountLSB,                0u, bitCountMSB);
}

void leb__SetData(uint nodeID, int nodeDepth, uint bitData)
{
    int bitCount = leb__DataBitSize(nodeDepth);

    leb__SetDataExplicit(nodeID, nodeDepth, bitCount, bitData);
}


/*******************************************************************************
 * GetData -- Returns bitCount bits located at arrayID
 *
 * Note that this procedure writes to at most two uint32 elements.
 * Two elements are relevant whenever the specified interval overflows 32-bit
 * words.
 *
 */
uint leb__GetDataExplicit(uint nodeID, int nodeDepth, uint bitCount)
{
    uint alignedBitOffset = leb__DataBitID(nodeID, nodeDepth);
    uint maxBufferIndex = leb__BufferUintSize(LEB_MAX_DEPTH) - 1;
    uint lebBufferIndexLSB = (alignedBitOffset >> 5u);
    uint lebBufferIndexMSB = min(lebBufferIndexLSB + 1, maxBufferIndex);
    uint bitFieldOffsetLSB = alignedBitOffset & 31u;
    uint bitCountLSB = min(32u - bitFieldOffsetLSB, bitCount);
    uint bitCountMSB = bitCount - bitCountLSB;
    uint bitFieldLSB = LEB_BUFFER[lebBufferIndexLSB];
    uint bitFieldMSB = LEB_BUFFER[lebBufferIndexMSB];
    uint lsb = leb__BitFieldExtract(bitFieldLSB, bitFieldOffsetLSB, bitCountLSB);
    uint msb = leb__BitFieldExtract(bitFieldMSB,                0u, bitCountMSB);

    return (lsb | (msb << bitCountLSB));
}

uint leb__GetData(uint nodeID, int nodeDepth)
{
    int bitCount = leb__DataBitSize(nodeDepth);

    return leb__GetDataExplicit(nodeID, nodeDepth, bitCount);
}


/*******************************************************************************
 * IsLeafNode -- Checks if a node is a leaf node
 *
 */
bool leb_IsLeafNode(in const leb_Node n)
{
    return (n.depth == LEB_MAX_DEPTH);
}


/*******************************************************************************
 * IsRootNode -- Checks if a node is a root node
 *
 */
bool leb_IsRootNode(in const leb_Node n)
{
    return (n.depth == 1);
}


/*******************************************************************************
 * IsNullNode -- Checks if a node is a null node
 *
 */
bool leb_IsNullNode(in const leb_Node n)
{
    return (n.id == 0u);
}


/*******************************************************************************
 * ParentNode -- Computes the parent of the input node
 *
 */
leb_Node leb_ParentNode(in const leb_Node node)
{
    return leb_Node(node.id >> 1u, node.depth - 1);
}


/*******************************************************************************
 * SiblingNode -- Computes the sibling of the input node
 *
 */
leb_Node leb__SiblingNode(in const leb_Node node)
{
    if (leb_IsNullNode(node)) {
        return leb_Node(0u, 0);
    } else {
        return leb_Node(node.id ^ 1u, node.depth);
    }
}


/*******************************************************************************
 * RightSiblingNode -- Computes the right sibling of the input node
 *
 */
leb_Node leb__RightSiblingNode(in const leb_Node node)
{
    if (leb_IsNullNode(node)) {
        return leb_Node(0u, 0);
    } else {
        return leb_Node(node.id | 1u, node.depth);
    }
}


/*******************************************************************************
 * LeftSiblingNode -- Computes the left sibling of the input node
 *
 */
leb_Node leb__LeftSiblingNode(in const leb_Node node)
{
    if (leb_IsNullNode(node)) {
        return leb_Node(0u, 0);
    } else {
        return leb_Node(node.id & (~1u), node.depth);
    }
}


/*******************************************************************************
 * RightChildNode -- Computes the right child of the input node
 *
 */
leb_Node leb__RightChildNode(in const leb_Node node)
{
    if (leb_IsNullNode(node)) {
        return leb_Node(0u, 0);
    } else {
        return leb_Node(node.id << 1u | 1u, node.depth + 1);
    }
}


/*******************************************************************************
 * LeftChildNode -- Computes the left child of the input node
 *
 */
leb_Node leb__LeftChildNode(in const leb_Node node)
{
    if (leb_IsNullNode(node)) {
        return leb_Node(0u, 0);
    } else {
        return leb_Node(node.id << 1u, node.depth + 1);
    }
}


/*******************************************************************************
 * Split -- Subdivides a node in two
 *
 */
void leb__SplitNode(in const leb_Node node)
{
    leb__SetNodeBitValue(leb__RightChildNode(node), 1u);
}

/*******************************************************************************
 * Merge -- Merges the node with its neighbour
 *
 */
void leb__MergeNode(in const leb_Node node)
{
    leb__SetNodeBitValue(leb__RightSiblingNode(node), 0u);
}


/*******************************************************************************
 * NodeCount -- Returns the number of triangles in the LEB
 *
 */
uint leb_NodeCount()
{
    return leb__GetData(1u, 0);
}


/*******************************************************************************
 * Decode the LEB Node associated to an index
 *
 */
leb_Node leb_DecodeNode(uint nodeID)
{
    leb_Node node = {1u, 0};

    while (leb__GetData(node.id, node.depth) > 1u) {
        uint cmp = leb__GetData(node.id<<= 1u, ++node.depth);
        uint b = nodeID < cmp ? 0 : 1;

        node.id|= b;
        nodeID-= cmp * b;
    }

    return node;
}


/*******************************************************************************
 * EncodeNode -- Returns the bit index associated with the Node
 *
 * This does the inverse of the DecodeNode routine.
 *
 */
uint leb_EncodeNode(in const leb_Node node)
{
    uint nodeID = 0u;
    leb_Node nodeIterator = node;

    while (nodeIterator.id > 1u) {
        leb_Node sibling = leb__LeftSiblingNode(nodeIterator);
        uint nodeCount = leb__GetData(sibling.id, sibling.depth);

        nodeID+= (nodeIterator.id & 1u) * nodeCount;
        nodeIterator = leb_ParentNode(nodeIterator);
    }

    return nodeID;
}


/*******************************************************************************
 * SplitNodeIDs -- Updates the IDs of neighbors after one LEB split
 *
 * This code applies the following rules:
 * Split left:
 * LeftID  = 2 * NodeID + 1
 * RightID = 2 * EdgeID + 1
 * EdgeID  = 2 * RightID + 1
 *
 * Split right:
 * LeftID  = 2 * EdgeID
 * RightID = 2 * NodeID
 * EdgeID  = 2 * LeftID
 *
 * The _reserved channel stores NodeID, which is recquired for applying the
 * rules.
 *
 */
leb_SameDepthNeighborIDs
leb__SplitNodeIDs(in const leb_SameDepthNeighborIDs nodeIDs, uint splitBit)
{
#if 1 // branchless version
    uint b = splitBit;
    uint c = splitBit ^ 1u;
    bool cb = bool(c);
    uvec4 idArray = uvec4(nodeIDs.left, nodeIDs.right, nodeIDs.edge, nodeIDs._reserved);
    leb_SameDepthNeighborIDs newIDs = {
        (idArray[2 + b] << 1u) | uint(cb && bool(idArray[2 + b])),
        (idArray[2 + c] << 1u) | uint(cb && bool(idArray[2 + c])),
        (idArray[b    ] << 1u) | uint(cb && bool(idArray[b    ])),
        (idArray[3    ] << 1u) | b
    };

    return newIDs;
#else
    uint n1 = nodeIDs.left, n2 = nodeIDs.right,
         n3 = nodeIDs.edge, n4 = nodeIDs._reserved;
    uint b2 = (n2 == 0u) ? 0u : 1u,
         b3 = (n3 == 0u) ? 0u : 1u;

    if (splitBit == 0u) {
        return leb_SameDepthNeighborIDs(
            n4 << 1 | 1, n3 << 1 | b3, n2 << 1 | b2, n4 << 1
        );
    } else {
        return leb_SameDepthNeighborIDs(
            n3 << 1    , n4 << 1     , n1 << 1     , n4 << 1 | 1
        );
    }
#endif
}

/*******************************************************************************
 * DecodeNodeNeighborIDs -- Decodes the IDs of the leb_Nodes neighbour to node
 *
 * The IDs are associated to the depth of the input node. As such, they
 * don't necessarily exist in the LEB subdivision.
 *
 */
leb_SameDepthNeighborIDs leb_DecodeSameDepthNeighborIDs(in const leb_Node node)
{
    uint b = leb__GetBitValue(node.id, node.depth - 1);
    leb_SameDepthNeighborIDs nodeIDs = leb_SameDepthNeighborIDs(0u, 0u, 3u - b, 2u + b);

    for (int bitID = node.depth - 2; bitID >= 0; --bitID) {
        nodeIDs = leb__SplitNodeIDs(nodeIDs, leb__GetBitValue(node.id, bitID));
    }

    return nodeIDs;
}


/*******************************************************************************
 * SameDepthNeighborIDs -- Computes the IDs of the same-level neighbors of a node
 *
 */
leb_SameDepthNeighborIDs
leb_GetSameDepthNeighborIDs(in const leb_NodeAndNeighbors nodes)
{
    uint edgeID = nodes.edge.id << (nodes.node.depth - nodes.edge.depth);
    uint leftID = nodes.left.id >> (nodes.left.depth - nodes.node.depth);
    uint rightID = nodes.right.id >> (nodes.right.depth - nodes.node.depth);

    return leb_SameDepthNeighborIDs(leftID, rightID, edgeID, nodes.node.id);
}


/*******************************************************************************
 * EdgeNode -- Computes the neighbour of the input node wrt to its longest edge
 *
 */
leb_Node EdgeNode(in const leb_Node node)
{
    return leb_Node(leb_DecodeSameDepthNeighborIDs(node).edge, node.depth);
}


/*******************************************************************************
 * SplitNodeConforming -- Splits a node while producing a conforming LEB
 *
 */
void leb_SplitNodeConforming(in const leb_Node node)
{
    if (!leb_IsLeafNode(node)) {
        leb_Node nodeIterator = node;

        leb__SplitNode(nodeIterator);
        nodeIterator = EdgeNode(nodeIterator);

        while (nodeIterator.id > 3u) {
            leb__SplitNode(nodeIterator);
            nodeIterator = leb_ParentNode(nodeIterator);
            leb__SplitNode(nodeIterator);
            nodeIterator = EdgeNode(nodeIterator);
        }
    }
}


/*******************************************************************************
 * HasNode -- Checks if the specified node really exists in the LEB
 *
 */
bool leb__HasNode(in const leb_Node node)
{
    uint nodeBit = leb__GetNodeBitValue(node)
                 & leb__GetNodeBitValue(leb__SiblingNode(node));

    if (!leb_IsLeafNode(node)) {
        nodeBit&= 1u ^ leb__GetNodeBitValue(leb__RightChildNode(node));
    }

    return (nodeBit == 1u);
}


/*******************************************************************************
 * MergeNodeConforming -- Merges a node while producing a conforming LEB
 *
 * This routines makes sure that the children of a diamond (including the
 * input node) all exist in the LEB before calling a merge.
 *
 */
void
leb_MergeNodeConforming(in const leb_Node node, in const leb_DiamondParent diamond)
{
    if (!leb_IsRootNode(node)) {
        uint id1 = node.id ^ 1u;
        uint id2 = diamond.top.id << 1u;
        uint id3 = diamond.top.id << 1u | 1u;
        bool b1 = leb__HasNode(leb_Node(id1, node.depth));
        bool b2 = leb__HasNode(leb_Node(id2, node.depth));
        bool b3 = leb__HasNode(leb_Node(id3, node.depth));

        if (b1 && b2 && b3) {
            leb__MergeNode(node);
            leb__MergeNode(leb_Node(id3, node.depth));
        }
    }
}


/*******************************************************************************
 * DecodeNodeDiamondIDs -- Decodes the upper Diamond associated to the leb_Node
 *
 * If the neighbour part does not exist, the parentNode is copied instead.
 *
 */
leb_DiamondParent leb_DecodeDiamondParent(in const leb_Node node)
{
    leb_Node parentNode = leb_ParentNode(node);
    uint diamondNodeID = leb_DecodeSameDepthNeighborIDs(parentNode).edge;
    leb_Node diamondNode = leb_Node(
        diamondNodeID > 0u ? diamondNodeID : parentNode.id,
        parentNode.depth
    );

    return leb_DiamondParent(parentNode, diamondNode);
}


/*******************************************************************************
 * SplitMatrix3x3 -- Computes a LEB splitting matrix from a split bit
 *
 */
mat3 leb__SplitMatrix3x3(uint splitBit)
{
    float b = float(splitBit);
    float c = 1.0f - b;

    return transpose(mat3(
        c   , b   , 0.0f,
        0.5f, 0.0f, 0.5f,
        0.0f,    c,    b
    ));
}


/*******************************************************************************
 * QuadMatrix3x3 -- Computes the matrix that affects the triangle to the quad
 *
 */
mat3 leb__QuadMatrix3x3(uint quadBit)
{
    float b = float(quadBit);
    float c = 1.0f - b;

    return transpose(mat3(
        c, 0.0f, b,
        b, c   , b,
        b, 0.0f, c
    ));
}


/*******************************************************************************
 * DecodeSplittingMatrix -- Computes the splitting matrix associated to a LEB
 * node
 *
 */
mat3 leb__DecodeSplittingMatrix(in const leb_Node node)
{
    int bitID = node.depth - 1;
    mat3 matrix = leb__QuadMatrix3x3(leb__GetBitValue(node.id, bitID));

    --bitID;

    while (bitID >= 0) {
        matrix = leb__SplitMatrix3x3(leb__GetBitValue(node.id, bitID)) * matrix;
        --bitID;
    }

    return matrix;
}


/*******************************************************************************
 * DecodeNodeAttributeArray -- Compute the triangle attributes at the input node
 *
 */
vec3 leb_DecodeNodeAttributeArray(in const leb_Node node, in const vec3 data)
{
    return leb__DecodeSplittingMatrix(node) * data;
}

mat2x3 leb_DecodeNodeAttributeArray(in const leb_Node node, in const mat2x3 data)
{
    return leb__DecodeSplittingMatrix(node) * data;
}

mat3x3 leb_DecodeNodeAttributeArray(in const leb_Node node, in const mat3x3 data)
{
    return leb__DecodeSplittingMatrix(node) * data;
}

mat4x3 leb_DecodeNodeAttributeArray(in const leb_Node node, in const mat4x3 data)
{
    return leb__DecodeSplittingMatrix(node) * data;
}


/*******************************************************************************
 * BoundingNode -- Compute the triangle that bounds the point (x, y)
 *
 */
leb_Node leb_BoundingNode(vec2 p)
{
    leb_Node node = leb_Node(0u, 0);

    if (p.x >= 0.0f && p.y >= 0.0f && p.x + p.y <= 1.0f) {
        node = leb_Node(1u, 0);

        while (!leb__HasNode(node)) {
            vec2 q = p;

            if (q.x < q.y) {
                node = leb__LeftChildNode(node);
                p.x = (1.0f - q.x - q.y);
                p.y = (q.y - q.x);
            } else {
                node = leb__RightChildNode(node);
                p.x = (q.x - q.y);
                p.y = (1.0f - q.x - q.y);
            }
        }
    }

    return node;
}

