/* leb.h - public domain Longest Edge Bisection library
by Jonathan Dupuy

   This is a low-level library for computing longest-edge bisections.

   INTERFACING

   define LEB_ASSERT(x) to avoid using assert.h.
   define LEB_LOG(format, ...) to use your own logger (default prints in stdout)
   define LEB_MALLOC(x) to use your own memory allocator
   define LEB_FREE(x) to use your own memory deallocator

*/

#ifndef LEB_INCLUDE_LEB_H
#define LEB_INCLUDE_LEB_H

#ifdef __cplusplus
extern "C" {
#endif

#ifdef LEB_STATIC
#define LEBDEF static
#else
#define LEBDEF extern
#endif

#include <stdint.h>

// data structures
typedef struct leb_Heap leb_Heap;
typedef struct {
    uint32_t id;
    int32_t depth;
} leb_Node;
typedef struct {
    uint32_t left, right, edge, _reserved;
} leb_SameDepthNeighborIDs;
typedef struct {
    leb_Node base, top;
} leb_DiamondParent;
typedef struct {
    leb_Node left, right, edge, node;
} leb_NodeAndNeighbors;

// ctor / dtor
// Note: maxDepth *must* be at least 5
LEBDEF leb_Heap *leb_Create(int maxDepth);
LEBDEF leb_Heap *leb_CreateMinMax(int minDepth, int maxDepth);
LEBDEF void leb_Release(leb_Heap *leb);

// serialization
LEBDEF const char  *leb_GetHeapMemory(const leb_Heap *leb);
LEBDEF void         leb_SetHeapMemory(leb_Heap *leb, const char *buffer);
LEBDEF uint32_t     leb_HeapByteSize(const leb_Heap *leb);

// loaders
LEBDEF void leb_ResetToRoot(leb_Heap *leb);
LEBDEF void leb_ResetToLeaf(leb_Heap *leb);
LEBDEF void leb_ResetToDepth(leb_Heap *leb, int depth);

// manipulation
LEBDEF void leb_ComputeSumReduction(leb_Heap *leb);
LEBDEF void leb_SplitNodeConforming(leb_Heap *leb, const leb_Node node);
LEBDEF void leb_MergeNodeConforming(leb_Heap *leb,
                                    const leb_Node node,
                                    const leb_DiamondParent diamond);
LEBDEF void leb_SplitNodeConforming_Quad(leb_Heap *leb, const leb_Node node);
LEBDEF void leb_MergeNodeConforming_Quad(leb_Heap *leb,
                                         const leb_Node node,
                                         const leb_DiamondParent diamond);

// O(1) queries
LEBDEF int32_t leb_MinDepth(const leb_Heap *leb);
LEBDEF int32_t leb_MaxDepth(const leb_Heap *leb);
LEBDEF uint32_t leb_NodeCount(const leb_Heap *leb);
LEBDEF bool leb_IsLeafNode(const leb_Heap *leb, const leb_Node node);
LEBDEF bool leb_IsCeilNode(const leb_Heap *leb, const leb_Node node);
LEBDEF bool leb_IsRootNode(const leb_Heap *leb, const leb_Node node);
LEBDEF bool leb_IsNullNode(                     const leb_Node node);
LEBDEF leb_Node leb_ParentNode(const leb_Node node);
LEBDEF leb_SameDepthNeighborIDs
       leb_GetSameDepthNeighborIDs(const leb_NodeAndNeighbors nodes);

// O(depth) queries
LEBDEF leb_Node leb_DecodeNode(const leb_Heap *leb, uint32_t bitID);
LEBDEF uint32_t leb_EncodeNode(const leb_Heap *leb, const leb_Node node);
LEBDEF leb_NodeAndNeighbors leb_DecodeNodeAndNeighbors(const leb_Heap *leb,
                                                       uint32_t bitID);
LEBDEF leb_NodeAndNeighbors leb_DecodeNodeAndNeighbors_Quad(const leb_Heap *leb,
                                                            uint32_t bitID);
LEBDEF leb_SameDepthNeighborIDs leb_DecodeSameDepthNeighborIDs(const leb_Node node);
LEBDEF leb_SameDepthNeighborIDs leb_DecodeSameDepthNeighborIDs_Quad(const leb_Node node);
LEBDEF leb_DiamondParent leb_DecodeDiamondParent(const leb_Node node);
LEBDEF leb_DiamondParent leb_DecodeDiamondParent_Quad(const leb_Node node);

// subdivision routine O(depth)
LEBDEF void leb_DecodeNodeAttributeArray(const leb_Node node,
                                         int attributeArraySize,
                                         float attributeArray[][3]);
LEBDEF void leb_DecodeNodeAttributeArray_Quad(const leb_Node node,
                                              int attributeArraySize,
                                              float attributeArray[][3]);

// intersection test O(depth)
LEBDEF leb_Node leb_BoundingNode     (const leb_Heap *leb, float x, float y);
LEBDEF leb_Node leb_BoundingNode_Quad(const leb_Heap *leb, float x, float y);


#ifdef __cplusplus
} // extern "C"
#endif

//
//
//// end header file ///////////////////////////////////////////////////////////
#endif // LEB_INCLUDE_LEB_H

#ifdef LEB_IMPLEMENTATION

#ifndef LEB_ASSERT
#    include <assert.h>
#    define LEB_ASSERT(x) assert(x)
#endif

#ifndef LEB_LOG
#    include <stdio.h>
#    define LEB_LOG(format, ...) do { fprintf(stdout, format, ##__VA_ARGS__); fflush(stdout); } while(0)
#endif

#ifndef LEB_MALLOC
#    include <stdlib.h>
#    define LEB_MALLOC(x) (malloc(x))
#    define LEB_FREE(x) (free(x))
#else
#    ifndef LEB_FREE
#        error LEB_MALLOC defined without LEB_FREE
#    endif
#endif


/*******************************************************************************
 * MinValue -- Returns the minimum value between two inputs
 *
 */
static inline uint32_t leb__MinValue(uint32_t a, uint32_t b)
{
    return a < b ? a : b;
}


/*******************************************************************************
 * GetBitValue -- Returns the value of a bit stored in a 32-bit word
 *
 */
static uint32_t leb__GetBitValue(const uint32_t bitField, uint32_t bitID)
{
    return ((bitField >> bitID) & 1u);
}


/*******************************************************************************
 * SetBitValue -- Sets the value of a bit stored in a 32-bit word
 *
 */
static void
leb__SetBitValue(uint32_t *bitField, uint32_t bitID, uint32_t bitValue)
{
    const uint32_t bitMask = ~(1u << bitID);

    (*bitField) = (*bitField & bitMask) | (bitValue << bitID);
}


/*******************************************************************************
 * BitFieldInsert -- Returns a uint32 bit field after insertion of some bit
 * data in range [bitOffset, bitOffset + bitCount - 1]
 *
 */
static inline void
leb__BitFieldInsert(
    uint32_t *bitField,
    uint32_t bitOffset,
    uint32_t bitCount,
    uint32_t bitData
) {
    LEB_ASSERT(bitOffset < 32u && bitCount <= 32u && bitOffset + bitCount <= 32u);
    uint32_t bitMask = ~(~(0xFFFFFFFFu << bitCount) << bitOffset);

    (*bitField) = (*bitField & bitMask) | (bitData << bitOffset);
}


/*******************************************************************************
 * BitFieldExtract -- Extracts bits [bitOffset, bitOffset + bitCount - 1] from
 * a uint32 bit field, returning them in the least significant bits of the result.
 *
 */
static inline uint32_t
leb__BitFieldExtract(
    const uint32_t bitField,
    uint32_t bitOffset,
    uint32_t bitCount
) {
    LEB_ASSERT(bitOffset < 32u && bitCount < 32u && bitOffset + bitCount <= 32u);
    uint32_t bitMask = ~(0xFFFFFFFFu << bitCount);

    return (bitField >> bitOffset) & bitMask;
}


/*******************************************************************************
 * Leb Buffer Data structure
 *
 */
struct leb_Heap {
    uint32_t *buffer;
    int32_t minDepth, maxDepth;
};


/*******************************************************************************
 * IsCeilNode -- Checks if a node is a ceil node, i.e., that can not split further
 *
 */
LEBDEF bool leb_IsCeilNode(const leb_Heap *leb, const leb_Node n)
{
    return (n.depth == leb->maxDepth);
}


/*******************************************************************************
 * IsRootNode -- Checks if a node is a root node
 *
 */
LEBDEF bool leb_IsRootNode(const leb_Heap *leb, const leb_Node n)
{
    return (n.depth == leb->minDepth);
}


/*******************************************************************************
 * IsNullNode -- Checks if a node is a null node
 *
 */
LEBDEF bool leb_IsNullNode(const leb_Node n)
{
    return (n.id == (uint32_t)(n.depth) /* == 0*/);
}


/*******************************************************************************
 * ParentNode -- Computes the parent of the input node
 *
 */
static leb_Node leb__ParentNode_Fast(const leb_Node node)
{
    return {node.id >> 1u, node.depth - 1};
}

LEBDEF leb_Node leb_ParentNode(const leb_Node node)
{
     return leb_IsNullNode(node) ? node : leb__ParentNode_Fast(node);
}


/*******************************************************************************
 * CeilNode -- Returns the associated ceil node, i.e., the deepest possible leaf
 *
 */
static leb_Node leb__CeilNode_Fast(const leb_Heap *leb, const leb_Node node)
{
    return {node.id << (leb->maxDepth - node.depth), leb->maxDepth};
}

static leb_Node leb__CeilNode(const leb_Heap *leb, const leb_Node node)
{
    return leb_IsNullNode(node) ? node : leb__CeilNode_Fast(leb, node);
}


/*******************************************************************************
 * SiblingNode -- Computes the sibling of the input node
 *
 */
static leb_Node leb__SiblingNode_Fast(const leb_Node node)
{
    return {node.id ^ 1u, node.depth};
}

static leb_Node leb__SiblingNode(const leb_Node node)
{
    return leb_IsNullNode(node) ? node : leb__SiblingNode_Fast(node);
}


/*******************************************************************************
 * RightSiblingNode -- Computes the right sibling of the input node
 *
 */
static leb_Node leb__RightSiblingNode_Fast(const leb_Node node)
{
    return {node.id | 1u, node.depth};
}

static leb_Node leb__RightSiblingNode(const leb_Node node)
{
    return leb_IsNullNode(node) ? node : leb__RightSiblingNode_Fast(node);
}


/*******************************************************************************
 * LeftSiblingNode -- Computes the left sibling of the input node
 *
 */
static leb_Node leb__LeftSiblingNode_Fast(const leb_Node node)
{
    return {node.id & (~1u), node.depth};
}

static leb_Node leb__LeftSiblingNode(const leb_Node node)
{
    return leb_IsNullNode(node) ? node : leb__LeftSiblingNode_Fast(node);
}


/*******************************************************************************
 * RightChildNode -- Computes the right child of the input node
 *
 */
static leb_Node leb__RightChildNode_Fast(const leb_Node node)
{
    return {node.id << 1u | 1u, node.depth + 1};
}

static leb_Node leb__RightChildNode(const leb_Node node)
{
    return leb_IsNullNode(node) ? node : leb__RightChildNode_Fast(node);
}


/*******************************************************************************
 * LeftChildNode -- Computes the left child of the input node
 *
 */
static leb_Node leb__LeftChildNode_Fast(const leb_Node node)
{
    return {node.id << 1u, node.depth + 1};
}

static leb_Node leb__LeftChildNode(const leb_Node node)
{
    return leb_IsNullNode(node) ? node : leb__LeftChildNode_Fast(node);
}


/*******************************************************************************
 * HeapBitSize -- Computes the number of bits to allocate for the buffer
 *
 * For a tree of max depth D, the number of bits is 2^(D+2).
 * Note that 2 bits are "wasted" in the sense that they only serve
 * to round the required number of bits to a power of two.
 *
 */
static inline uint32_t leb__HeapBitSize(uint32_t lebMaxDepth)
{
    return 1u << (lebMaxDepth + 2u);
}


/*******************************************************************************
 * HeapUint32Size -- Computes the number of uints to allocate for the bitfield
 *
 */
static inline uint32_t leb__HeapUint32Size(uint32_t lebMaxDepth)
{
    return leb__HeapBitSize(lebMaxDepth) >> 5u;
}


/*******************************************************************************
 * HeapByteSize -- Computes the number of Bytes to allocate for the bitfield
 *
 */
static uint32_t leb__HeapByteSize(uint32_t lebMaxDepth)
{
    return leb__HeapUint32Size(lebMaxDepth) * sizeof(uint32_t);
}


/*******************************************************************************
 * NodeBitID -- Returns the bit index that stores data associated with a given node
 *
 * For a LEB of max depth D and given an index in [0, 2^(D+1) - 1], this
 * functions is used to emulate the behaviour of a lookup in an array, i.e.,
 * uint32_t[nodeID]. It provides the first bit in memory that stores
 * information associated with the element of index nodeID.
 *
 * For data located at level d, the bit offset is 2^d x (3 - d + D)
 * We then offset this quantity by the index by (nodeID - 2^d) x (D + 1 - d)
 * Note that the null index (nodeID = 0) is also supported.
 *
 */
static inline uint32_t leb__NodeBitID(const leb_Heap *leb, const leb_Node node)
{
    uint32_t tmp1 = 2u << node.depth;
    uint32_t tmp2 = 1u + (uint32_t)(leb->maxDepth - node.depth);

    return tmp1 + node.id * tmp2;
}


/*******************************************************************************
 * NodeBitID_BitField -- Computes the bitfield bit location associated to a node
 *
 * Here, the node is converted into a final node and its bit offset is
 * returned, which is finalNodeID + 2^{D + 1}
 */
static uint32_t
leb__NodeBitID_BitField(const leb_Heap *leb, const leb_Node node)
{
    return leb__NodeBitID(leb, leb__CeilNode(leb, node));
}


/*******************************************************************************
 * NodeBitSize -- Returns the number of bits storing the input node value
 *
 */
static inline int32_t leb__NodeBitSize(const leb_Heap *leb, const leb_Node node)
{
    return leb->maxDepth - node.depth + 1;
}


/*******************************************************************************
 * HeapArgs
 *
 * The LEB heap data structure uses an array of 32-bit words to store its data.
 * Whenever we need to access a certain bit range, we need to query two such
 * words (because sometimes the requested bit range overlaps two 32-bit words).
 * The HeapArg data structure provides arguments for reading from and/or
 * writing to the two 32-bit words that bound the queries range.
 *
 */
typedef struct {
    uint32_t *bitFieldLSB, *bitFieldMSB;
    uint32_t bitOffsetLSB;
    uint32_t bitCountLSB, bitCountMSB;
} leb__HeapArgs;

leb__HeapArgs
leb__CreateHeapArgs(const leb_Heap *leb, const leb_Node node, int32_t bitCount)
{
    uint32_t alignedBitOffset = leb__NodeBitID(leb, node);
    uint32_t maxBufferIndex = leb__HeapUint32Size(leb->maxDepth) - 1u;
    uint32_t lebBufferIndexLSB = (alignedBitOffset >> 5u);
    uint32_t lebBufferIndexMSB = leb__MinValue(lebBufferIndexLSB + 1, maxBufferIndex);
    leb__HeapArgs args;

    args.bitOffsetLSB = alignedBitOffset & 31u;
    args.bitCountLSB = leb__MinValue(32u - args.bitOffsetLSB, bitCount);
    args.bitCountMSB = bitCount - args.bitCountLSB;
    args.bitFieldLSB = &leb->buffer[lebBufferIndexLSB];
    args.bitFieldMSB = &leb->buffer[lebBufferIndexMSB];

    return args;
}


/*******************************************************************************
 * HeapWrite -- Sets bitCount bits located at nodeID to bitData
 *
 * Note that this procedure writes to at most two uint32 elements.
 * Two elements are relevant whenever the specified interval overflows 32-bit
 * words.
 *
 */
static void
leb__HeapWriteExplicit(
    leb_Heap *leb,
    const leb_Node node,
    int32_t bitCount,
    uint32_t bitData
) {
    leb__HeapArgs args = leb__CreateHeapArgs(leb, node, bitCount);

    leb__BitFieldInsert(args.bitFieldLSB,
                        args.bitOffsetLSB,
                        args.bitCountLSB,
                        bitData);
    leb__BitFieldInsert(args.bitFieldMSB,
                        0u,
                        args.bitCountMSB,
                        bitData >> args.bitCountLSB);
}

static void leb__HeapWrite(leb_Heap *leb, const leb_Node node, uint32_t bitData)
{
    leb__HeapWriteExplicit(leb, node, leb__NodeBitSize(leb, node), bitData);
}


/*******************************************************************************
 * HeapRead -- Returns bitCount bits located at nodeID
 *
 * Note that this procedure writes to at most two uint32 elements.
 * Two elements are relevant whenever the specified interval overflows 32-bit
 * words.
 *
 */
static uint32_t
leb__HeapReadExplicit(
    const leb_Heap *leb,
    const leb_Node node,
    int32_t bitCount
) {
    leb__HeapArgs args = leb__CreateHeapArgs(leb, node, bitCount);
    uint32_t lsb = leb__BitFieldExtract(*args.bitFieldLSB,
                                        args.bitOffsetLSB,
                                        args.bitCountLSB);
    uint32_t msb = leb__BitFieldExtract(*args.bitFieldMSB,
                                        0u,
                                        args.bitCountMSB);

    return (lsb | (msb << args.bitCountLSB));
}

static uint32_t leb__HeapRead(const leb_Heap *leb, const leb_Node node)
{
    return leb__HeapReadExplicit(leb, node, leb__NodeBitSize(leb, node));
}


/*******************************************************************************
 * HeapWrite_BitField -- Sets the bit associated to a leaf node to bitValue
 *
 * This is a dedicated routine to write directly to the bitfield.
 *
 */
static void
leb__HeapWrite_BitField(
    leb_Heap *leb,
    const leb_Node node,
    const uint32_t bitValue
) {
    uint32_t bitID = leb__NodeBitID_BitField(leb, node);

    leb__SetBitValue(&leb->buffer[bitID >> 5u], bitID & 31u, bitValue);
}


/*******************************************************************************
 * HeapRead_BitField -- Returns the value of the bit associated to a leaf node
 *
 * This is a dedicated routine to read directly from the bitfield.
 *
 */
static uint32_t leb__HeapRead_BitField(const leb_Heap *leb, const leb_Node node)
{
    uint32_t bitID = leb__NodeBitID_BitField(leb, node);

    return leb__GetBitValue(leb->buffer[bitID >> 5u], bitID & 31u);
}


/*******************************************************************************
 * IsLeafNode -- Checks if a node is a leaf node, i.e., that has no descendants
 *
 */
LEBDEF bool leb_IsLeafNode(const leb_Heap *leb, const leb_Node node)
{
    return (leb__HeapRead(leb, node) == 1u);
}


/*******************************************************************************
 * ClearData -- Resets the data buffer stored by a LEB
 *
 */
static void leb__ClearBuffer(leb_Heap *leb)
{
    memset(leb->buffer, 0, leb__HeapByteSize(leb->maxDepth));
}


/*******************************************************************************
 * GetHeapMemory -- Returns a read-only pointer to the heap memory
 *
 */
LEBDEF const char *leb_GetHeapMemory(const leb_Heap *leb)
{
    return (char *)leb->buffer;
}


/*******************************************************************************
 * SetHeapMemory -- Sets the heap memory from a read-only buffer
 *
 */
LEBDEF void leb_SetHeapMemory(leb_Heap *leb, const char *buffer)
{
    memcpy(leb->buffer, buffer, leb_HeapByteSize(leb));
}


/*******************************************************************************
 * HeapByteSize -- Returns the amount of bytes consumed by the LEB heap
 *
 */
LEBDEF uint32_t leb_HeapByteSize(const leb_Heap *leb)
{
    return leb__HeapByteSize(leb->maxDepth);
}


/*******************************************************************************
 * Buffer Ctor
 *
 */
LEBDEF leb_Heap *leb_CreateMinMax(int minDepth, int maxDepth)
{
    LEB_ASSERT(maxDepth >=  5 && "maxDepth must be at least 5");
    LEB_ASSERT(maxDepth <= 29 && "maxDepth must be at most 29");
    LEB_ASSERT(minDepth >=  0 && "minDepth must be at least 0");
    LEB_ASSERT(minDepth <= maxDepth && "minDepth must be less than maxDepth");
    leb_Heap *leb = (leb_Heap *)LEB_MALLOC(sizeof(*leb));

    leb->minDepth = minDepth;
    leb->maxDepth = maxDepth;
    leb->buffer = (uint32_t *)LEB_MALLOC(leb__HeapByteSize(maxDepth));
    LEB_ASSERT(leb->buffer != NULL && "Memory allocation failed");
    leb_ResetToRoot(leb);

    return leb;
}

LEBDEF leb_Heap *leb_Create(int maxDepth)
{
    return leb_CreateMinMax(0, maxDepth);
}


/*******************************************************************************
 * Buffer Dtor
 *
 */
LEBDEF void leb_Release(leb_Heap *leb)
{
    LEB_FREE(leb->buffer);
    LEB_FREE(leb);
}


/*******************************************************************************
 * ResetToDepth -- Initializes a LEB to its a specific subdivision level
 *
 */
LEBDEF void leb_ResetToDepth(leb_Heap *leb, int depth)
{
    LEB_ASSERT(depth >= leb->minDepth && "depth must be at least equal to minDepth");
    LEB_ASSERT(depth <= leb->maxDepth && "depth must be at most equal to maxDepth");
    uint32_t minNodeID = 1u << depth;
    uint32_t maxNodeID = 2u << depth;

    leb__ClearBuffer(leb);

    for (uint32_t nodeID = minNodeID; nodeID < maxNodeID; ++nodeID) {
        leb_Node node = {nodeID, depth};

        leb__HeapWrite_BitField(leb, node, 1u);
    }

    leb_ComputeSumReduction(leb);
}


/*******************************************************************************
 * ResetToRoot -- Initializes a LEB to its minimum subdivision level
 *
 */
LEBDEF void leb_ResetToRoot(leb_Heap *leb)
{
    leb_ResetToDepth(leb, leb->minDepth);
}


/*******************************************************************************
 * ResetToLeaf -- Initializes a LEB to its maximum subdivision level
 *
 */
LEBDEF void leb_ResetToLeaf(leb_Heap *leb)
{
    leb_ResetToDepth(leb, leb->maxDepth);
}


/*******************************************************************************
 * UpdateBuffer -- Sums the 2 elements below the current slot
 *
 */
LEBDEF void leb_ComputeSumReduction(leb_Heap *leb)
{
    int depth = leb->maxDepth;
    uint32_t minNodeID = (1u << depth);
    uint32_t maxNodeID = (2u << depth);

    // prepass: processes deepest levels in parallel
    for (uint32_t nodeID = minNodeID; nodeID < maxNodeID; nodeID+= 32u) {
        uint32_t alignedBitOffset = leb__NodeBitID(leb, {nodeID, depth});
        uint32_t bitField = leb->buffer[alignedBitOffset >> 5u];
        uint32_t bitData = 0u;

        // 2-bits
        bitField = (bitField & 0x55555555u) + ((bitField >> 1u) & 0x55555555u);
        bitData = bitField;
        leb->buffer[(alignedBitOffset - minNodeID) >> 5u] = bitData;

        // 3-bits
        bitField = (bitField & 0x33333333u) + ((bitField >>  2u) & 0x33333333u);
        bitData = ((bitField >> 0u) & (7u <<  0u))
                | ((bitField >> 1u) & (7u <<  3u))
                | ((bitField >> 2u) & (7u <<  6u))
                | ((bitField >> 3u) & (7u <<  9u))
                | ((bitField >> 4u) & (7u << 12u))
                | ((bitField >> 5u) & (7u << 15u))
                | ((bitField >> 6u) & (7u << 18u))
                | ((bitField >> 7u) & (7u << 21u));
        leb__HeapWriteExplicit(leb, {nodeID >> 2u, depth - 2}, 24u, bitData);

        // 4-bits
        bitField = (bitField & 0x0F0F0F0Fu) + ((bitField >>  4u) & 0x0F0F0F0Fu);
        bitData = ((bitField >>  0u) & (15u <<  0u))
                | ((bitField >>  4u) & (15u <<  4u))
                | ((bitField >>  8u) & (15u <<  8u))
                | ((bitField >> 12u) & (15u << 12u));
        leb__HeapWriteExplicit(leb, {nodeID >> 3u, depth - 3}, 16u, bitData);

        // 5-bits
        bitField = (bitField & 0x00FF00FFu) + ((bitField >>  8u) & 0x00FF00FFu);
        bitData = ((bitField >>  0u) & (31u << 0u))
                | ((bitField >> 11u) & (31u << 5u));
        leb__HeapWriteExplicit(leb, {nodeID >> 4u, depth - 4}, 10u, bitData);

        // 6-bits
        bitField = (bitField & 0x0000FFFFu) + ((bitField >> 16u) & 0x0000FFFFu);
        bitData = bitField;
        leb__HeapWriteExplicit(leb, {nodeID >> 5u, depth - 5},  6u, bitData);
    }
    depth-= 5;

    // iterate over elements atomically
    while (--depth >= 0) {
        uint32_t minNodeID = 1u << depth;
        uint32_t maxNodeID = 2u << depth;

        for (uint32_t j = minNodeID; j < maxNodeID; ++j) {
            uint32_t x0 = leb__HeapRead(leb, {j << 1u     , depth + 1});
            uint32_t x1 = leb__HeapRead(leb, {j << 1u | 1u, depth + 1});

            leb__HeapWrite(leb, {j, depth}, x0 + x1);
        }
    }
}


/*******************************************************************************
 * Split -- Subdivides a node in two
 *
 */
static void leb__SplitNode(leb_Heap *leb, const leb_Node node)
{
    leb__HeapWrite_BitField(leb, leb__RightChildNode(node), 1u);
}


/*******************************************************************************
 * Merge -- Merges the node with its neighbour
 *
 */
static void leb__MergeNode(leb_Heap *leb, const leb_Node node)
{
    leb__HeapWrite_BitField(leb, leb__RightSiblingNode(node), 0u);
}


/*******************************************************************************
 * MinDepth -- Returns the min LEB depth
 *
 */
LEBDEF int32_t leb_MinDepth(const leb_Heap *leb)
{
    return leb->minDepth;
}


/*******************************************************************************
 * MaxDepth -- Returns the max LEB depth
 *
 */
LEBDEF int32_t leb_MaxDepth(const leb_Heap *leb)
{
    return leb->maxDepth;
}


/*******************************************************************************
 * NodeCount -- Returns the number of triangles in the LEB
 *
 */
LEBDEF uint32_t leb_NodeCount(const leb_Heap *leb)
{
    return leb__HeapRead(leb, {1u, 0});
}


/*******************************************************************************
 * DecodeNode -- Returns the LEB Node associated to index nodeID
 *
 * This is procedure is for iterating over the nodes.
 *
 */
LEBDEF leb_Node leb_DecodeNode(const leb_Heap *leb, uint32_t nodeID)
{
    LEB_ASSERT(nodeID < leb_NodeCount(leb) && "nodeID > NodeCount(leb)");

    leb_Node node = {1u, 0};

    while (leb__HeapRead(leb, node) > 1u) {
        uint32_t cmp = leb__HeapRead(leb, {node.id<<= 1u, ++node.depth});
        uint32_t b = nodeID < cmp ? 0 : 1;

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
LEBDEF uint32_t leb_EncodeNode(const leb_Heap *leb, const leb_Node node)
{
    uint32_t nodeID = 0u;
    leb_Node nodeIterator = node;

    while (nodeIterator.id > 1u) {
        leb_Node sibling = leb__LeftSiblingNode(nodeIterator);
        uint32_t nodeCount = leb__HeapRead(leb, {sibling.id, sibling.depth});

        nodeID+= (nodeIterator.id & 1u) * nodeCount;
        nodeIterator = leb__ParentNode_Fast(nodeIterator);
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
 * The _reserved channel stores NodeID, which is required for applying the
 * rules.
 *
 */
static leb_SameDepthNeighborIDs
leb__SplitNodeIDs(
    const leb_SameDepthNeighborIDs nodeIDs,
    const uint32_t splitBit
) {
    uint32_t n1 = nodeIDs.left, n2 = nodeIDs.right,
             n3 = nodeIDs.edge, n4 = nodeIDs._reserved;
    uint32_t b2 = (n2 == 0u) ? 0u : 1u,
             b3 = (n3 == 0u) ? 0u : 1u;

    if (splitBit == 0u) {
        return {n4 << 1 | 1, n3 << 1 | b3, n2 << 1 | b2, n4 << 1    };
    } else {
        return {n3 << 1    , n4 << 1     , n1 << 1     , n4 << 1 | 1};
    }
}


/*******************************************************************************
 * DecodeSameDepthNeighborIDs -- Decodes the IDs of the leb_Nodes neighbour to node
 *
 * The IDs are associated to the depth of the input node. As such, they
 * don't necessarily exist in the LEB subdivision.
 */
LEBDEF leb_SameDepthNeighborIDs
leb_DecodeSameDepthNeighborIDs(const leb_Node node)
{
    leb_SameDepthNeighborIDs nodeIDs = {0u, 0u, 0u, 1u};

    for (int bitID = node.depth - 1; bitID >= 0; --bitID) {
        nodeIDs = leb__SplitNodeIDs(nodeIDs, leb__GetBitValue(node.id, bitID));
    }

    return nodeIDs;
}

LEBDEF leb_SameDepthNeighborIDs
leb_DecodeSameDepthNeighborIDs_Quad(const leb_Node node)
{
    if (node.depth == 0)
        return {0u, 0u, 0u, 1u};

    uint32_t b = leb__GetBitValue(node.id, node.depth - 1);
    leb_SameDepthNeighborIDs nodeIDs = {0u, 0u, 3u - b, 2u + b};

    for (int bitID = node.depth - 2; bitID >= 0; --bitID) {
        nodeIDs = leb__SplitNodeIDs(nodeIDs, leb__GetBitValue(node.id, bitID));
    }

    return nodeIDs;
}


/*******************************************************************************
 * SameDepthNeighborIDs -- Computes the IDs of the same-level neighbors of a node
 *
 */
LEBDEF leb_SameDepthNeighborIDs
leb_GetSameDepthNeighborIDs(const leb_NodeAndNeighbors nodes)
{
    uint32_t edgeID = nodes.edge.id << (nodes.node.depth - nodes.edge.depth);
    uint32_t leftID = nodes.left.id >> (nodes.left.depth - nodes.node.depth);
    uint32_t rightID = nodes.right.id >> (nodes.right.depth - nodes.node.depth);

    return {leftID, rightID, edgeID, nodes.node.id};
}


/*******************************************************************************
 * EdgeNode -- Computes the neighbour of the input node wrt to its longest edge
 *
 */
static leb_Node leb__EdgeNode(const leb_Node node)
{
    uint32_t nodeID = leb_DecodeSameDepthNeighborIDs(node).edge;

    return {nodeID, (nodeID == 0u) ? 0 : node.depth};
}

static leb_Node leb__EdgeNode_Quad(const leb_Node node)
{
    uint32_t nodeID = leb_DecodeSameDepthNeighborIDs_Quad(node).edge;

    return {nodeID, (nodeID == 0u) ? 0 : node.depth};
}


/*******************************************************************************
 * SplitNodeConforming -- Splits a node while producing a conforming LEB
 *
 */
LEBDEF void leb_SplitNodeConforming(leb_Heap *leb, const leb_Node node)
{
    if (!leb_IsCeilNode(leb, node)) {
        const uint32_t minNodeID = 1u << leb->minDepth;
        leb_Node nodeIterator = node;

        leb__SplitNode(leb, nodeIterator);
        nodeIterator = leb__EdgeNode(nodeIterator);

        while (nodeIterator.id >= minNodeID) {
            leb__SplitNode(leb, nodeIterator);
            nodeIterator = leb_ParentNode(nodeIterator);
            leb__SplitNode(leb, nodeIterator);
            nodeIterator = leb__EdgeNode(nodeIterator);
        }
    }
}

LEBDEF void leb_SplitNodeConforming_Quad(leb_Heap *leb, const leb_Node node)
{
    if (!leb_IsCeilNode(leb, node)) {
        const uint32_t minNodeID = 1u << leb->minDepth;
        leb_Node nodeIterator = node;

        leb__SplitNode(leb, nodeIterator);
        nodeIterator = leb__EdgeNode_Quad(nodeIterator);

        while (nodeIterator.id >= minNodeID) {
            leb__SplitNode(leb, nodeIterator);
            nodeIterator = leb_ParentNode(nodeIterator);
            leb__SplitNode(leb, nodeIterator);
            nodeIterator = leb__EdgeNode_Quad(nodeIterator);
        }
    }
}


/*******************************************************************************
 * MergeNodeConforming -- Merges a node while producing a conforming LEB
 *
 * This routines makes sure that the children of a diamond (including the
 * input node) all exist in the LEB before calling a merge.
 *
 */
LEBDEF void
leb_MergeNodeConforming(
    leb_Heap *leb,
    const leb_Node node,
    const leb_DiamondParent diamond
) {
    if (!leb_IsRootNode(leb, node)) {
        leb_Node dualNode = leb__RightChildNode(diamond.top);
        bool b1 = leb_IsLeafNode(leb, leb__SiblingNode_Fast(node));
        bool b2 = leb_IsLeafNode(leb, dualNode);
        bool b3 = leb_IsLeafNode(leb, leb__SiblingNode(dualNode));

        if (b1 && b2 && b3) {
            leb__MergeNode(leb, node);
            leb__MergeNode(leb, dualNode);
        }
    }
}

LEBDEF void
leb_MergeNodeConforming_Quad(
    leb_Heap *leb,
    const leb_Node node,
    const leb_DiamondParent diamond
) {
    leb_MergeNodeConforming(leb, node, diamond);
}


/*******************************************************************************
 * DecodeDiamondParent -- Decodes the upper Diamond associated to the leb_Node
 *
 * If the neighbour part does not exist, the parentNode is copied instead.
 *
 */
LEBDEF leb_DiamondParent leb_DecodeDiamondParent(const leb_Node node)
{
    leb_Node parentNode = leb_ParentNode(node);
    uint32_t diamondNodeID = leb_DecodeSameDepthNeighborIDs(parentNode).edge;
    leb_Node diamondNode = {
        diamondNodeID > 0u ? diamondNodeID : parentNode.id,
        parentNode.depth
    };

    return {parentNode, diamondNode};
}

LEBDEF leb_DiamondParent leb_DecodeDiamondParent_Quad(const leb_Node node)
{
    leb_Node parentNode = leb_ParentNode(node);
    uint32_t diamondNodeID = leb_DecodeSameDepthNeighborIDs_Quad(parentNode).edge;
    leb_Node diamondNode = {
        diamondNodeID > 0u ? diamondNodeID : parentNode.id,
        parentNode.depth
    };

    return {parentNode, diamondNode};
}


/*******************************************************************************
 * NodeAndNeighborsFromSameDepthNeighborIDs -- Decodes the true neighbors of a node
 *
 */
static leb_NodeAndNeighbors
leb__NodeAndNeighborsFromSameDepthNeighborIDs(
    const leb_Heap *leb,
    const leb_SameDepthNeighborIDs nodeIDs,
    int nodeDepth
) {
    leb_NodeAndNeighbors nodeData = {
        {nodeIDs.left, nodeDepth},
        {nodeIDs.right, nodeDepth},
        {nodeIDs.edge, nodeDepth},
        {nodeIDs._reserved, nodeDepth}
    };

    if (!leb_IsLeafNode(leb, nodeData.edge))
        nodeData.edge = leb_ParentNode(nodeData.edge);

    if (!leb_IsLeafNode(leb, nodeData.left))
        nodeData.left = leb__RightChildNode(nodeData.left);

    if (!leb_IsLeafNode(leb, nodeData.right))
        nodeData.right = leb__LeftChildNode(nodeData.right);

    return nodeData;
}


/*******************************************************************************
 * DecodeNodeAndNeighbors -- Decode the LEB Node associated to an index, along with its neighbors
 *
 */
LEBDEF leb_NodeAndNeighbors
leb_DecodeNodeAndNeighbors(const leb_Heap *leb, uint32_t threadID)
{
#define nodeID nodeIDs._reserved
    leb_SameDepthNeighborIDs nodeIDs = {0u, 0u, 0u, 1u};
    int32_t nodeDepth = 0;

    while (leb__HeapRead(leb, {nodeID, nodeDepth}) > 1u) {
        uint32_t cmp = leb__HeapRead(leb, {nodeID << 1u, ++nodeDepth});
        uint32_t b = threadID < cmp ? 0u : 1u;

        nodeIDs = leb__SplitNodeIDs(nodeIDs, b);
        threadID-= cmp * b;
    }

    return leb__NodeAndNeighborsFromSameDepthNeighborIDs(leb, nodeIDs, nodeDepth);
#undef nodeID
}


/******************************************************************************/
/* Standalone matrix 3x3 API
 *
 */
typedef float lebMatrix3x3[3][3];


/*******************************************************************************
 * IdentityMatrix3x3 -- Sets a 3x3 matrix to identity
 *
 */
static void leb__IdentityMatrix3x3(lebMatrix3x3 m)
{
    m[0][0] = 1.0f; m[0][1] = 0.0f; m[0][2] = 0.0f;
    m[1][0] = 0.0f; m[1][1] = 1.0f; m[1][2] = 0.0f;
    m[2][0] = 0.0f; m[2][1] = 0.0f; m[2][2] = 1.0f;
}


/*******************************************************************************
 * TransposeMatrix3x3 -- Transposes a 3x3 matrix
 *
 */
static void leb__TransposeMatrix3x3(const lebMatrix3x3 m, lebMatrix3x3 out)
{
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            out[i][j] = m[j][i];
}


/*******************************************************************************
 * DotProduct -- Returns the dot product of two vectors
 *
 */
static float leb__DotProduct(int argSize, const float *x, const float *y)
{
    float dp = 0.0f;

    for (int i = 0; i < argSize; ++i)
        dp+= x[i] * y[i];

    return dp;
}


/*******************************************************************************
 * MulMatrix3x3 -- Computes the product of two 3x3 matrices
 *
 */
static void
leb__Matrix3x3Product(
    const lebMatrix3x3 m1,
    const lebMatrix3x3 m2,
    lebMatrix3x3 out
) {
    lebMatrix3x3 tra;

    leb__TransposeMatrix3x3(m2, tra);

    for (int j = 0; j < 3; ++j)
    for (int i = 0; i < 3; ++i)
        out[j][i] = leb__DotProduct(3, m1[j], tra[i]);
}


/*******************************************************************************
 * SplitMatrix3x3 -- Computes a LEB splitting matrix from a split bit
 *
 */
static void
leb__SplittingMatrix(lebMatrix3x3 splittingMatrix, uint32_t splitBit)
{
    float b = float(splitBit);
    float c = 1.0f - b;
    lebMatrix3x3 splitMatrix = {
        {c   , b   , 0.0f},
        {0.5f, 0.0f, 0.5f},
        {0.0f,    c,    b}
    };
    lebMatrix3x3 tmp;

    memcpy(tmp, splittingMatrix, sizeof(tmp));
    leb__Matrix3x3Product(splitMatrix, tmp, splittingMatrix);
}


/*******************************************************************************
 * QuadMatrix3x3 -- Computes a mirroring matrix from a split bit
 *
 */
static void
leb__QuadMatrix(lebMatrix3x3 matrix, uint32_t splitBit)
{
    float b = float(splitBit);
    float c = 1.0f - b;
    lebMatrix3x3 quadMatrix = {
        {c, 0.0f,    b},
        {b,    c,    b},
        {b, 0.0f,    c}
    };

    memcpy(matrix, quadMatrix, sizeof(quadMatrix));
}


/*******************************************************************************
 * DecodeTransformationMatrix -- Computes the matrix associated to a LEB
 * node
 *
 */
static void
leb__DecodeTransformationMatrix(
    const leb_Node node,
    lebMatrix3x3 splittingMatrix
) {
    leb__IdentityMatrix3x3(splittingMatrix);

    for (int bitID = node.depth - 1; bitID >= 0; --bitID) {
        leb__SplittingMatrix(splittingMatrix, leb__GetBitValue(node.id, bitID));
    }
}

static void
leb__DecodeTransformationMatrix_Quad(
    const leb_Node node,
    lebMatrix3x3 splittingMatrix
) {
    leb__QuadMatrix(splittingMatrix, leb__GetBitValue(node.id, node.depth - 1));

    for (int bitID = node.depth - 2; bitID >= 0; --bitID) {
        leb__SplittingMatrix(splittingMatrix, leb__GetBitValue(node.id, bitID));
    }
}


/*******************************************************************************
 * DecodeNodeAttributeArray -- Compute the triangle attributes at the input node
 *
 */
LEBDEF void
leb_DecodeNodeAttributeArray(
    const leb_Node node,
    int attributeArraySize,
    float attributeArray[][3]
) {
    LEB_ASSERT(attributeArraySize > 0);

    lebMatrix3x3 m;
    float attributeVector[3];

    leb__DecodeTransformationMatrix(node, m);

    for (int i = 0; i < attributeArraySize; ++i) {
        memcpy(attributeVector, attributeArray[i], sizeof(attributeVector));
        attributeArray[i][0] = leb__DotProduct(3, m[0], attributeVector);
        attributeArray[i][1] = leb__DotProduct(3, m[1], attributeVector);
        attributeArray[i][2] = leb__DotProduct(3, m[2], attributeVector);
    }
}

LEBDEF void
leb_DecodeNodeAttributeArray_Quad(
    const leb_Node node,
    int attributeArraySize,
    float attributeArray[][3]
) {
    LEB_ASSERT(attributeArraySize > 0);

    lebMatrix3x3 m;
    float attributeVector[3];

    leb__DecodeTransformationMatrix_Quad(node, m);

    for (int i = 0; i < attributeArraySize; ++i) {
        memcpy(attributeVector, attributeArray[i], sizeof(attributeVector));
        attributeArray[i][0] = leb__DotProduct(3, m[0], attributeVector);
        attributeArray[i][1] = leb__DotProduct(3, m[1], attributeVector);
        attributeArray[i][2] = leb__DotProduct(3, m[2], attributeVector);
    }
}


/*******************************************************************************
 * BoundingNode -- Compute the triangle that bounds the point (x, y)
 *
 */
LEBDEF leb_Node leb_BoundingNode(const leb_Heap *leb, float x, float y)
{
    leb_Node node = {0u, 0};

    if (x >= 0.0f && y >= 0.0f && x + y <= 1.0f) {
        node = {1u, 0};

        while (!leb_IsLeafNode(leb, node) && !leb_IsCeilNode(leb, node)) {
            float s = x, t = y;

            if (s < t) {
                node = leb__LeftChildNode(node);
                x = (1.0f - s - t);
                y = (t - s);
            } else {
                node = leb__RightChildNode(node);
                x = (s - t);
                y = (1.0f - s - t);
            }
        }
    }

    return node;
}

LEBDEF leb_Node leb_BoundingNode_Quad(const leb_Heap *leb, float x, float y)
{
    leb_Node node = {0u, 0};

    if (x >= 0.0f && y >= 0.0f && x <= 1.0f && y <= 1.0f) {
        if (x + y <= 1.0f) {
            node = {2u, 1};
        } else {
            node = {3u, 1};
            x = 1 - x;
            y = 1 - y;
        }

        while (!leb_IsLeafNode(leb, node) && !leb_IsCeilNode(leb, node)) {
            float s = x, t = y;

            if (s < t) {
                node = leb__LeftChildNode(node);
                x = (1.0f - s - t);
                y = (t - s);
            } else {
                node = leb__RightChildNode(node);
                x = (s - t);
                y = (1.0f - s - t);
            }
        }
    }

    return node;
}

#endif // LEB_IMPLEMENTATION
