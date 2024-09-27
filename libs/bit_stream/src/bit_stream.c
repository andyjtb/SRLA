#include "bit_stream.h"
#include <stdint.h>

/* Mask to extract lower bits up to 32 bits */
const uint32_t g_bitstream_lower_bits_mask[33] = {
    0x00000000U,
    0x00000001U, 0x00000003U, 0x00000007U, 0x0000000FU,
    0x0000001FU, 0x0000003FU, 0x0000007FU, 0x000000FFU,
    0x000001FFU, 0x000003FFU, 0x000007FFU, 0x00000FFFU,
    0x00001FFFU, 0x00003FFFU, 0x00007FFFU, 0x0000FFFFU,
    0x0001FFFFU, 0x0003FFFFU, 0x0007FFFFU, 0x000FFFFFU,
    0x001FFFFFU, 0x003FFFFFU, 0x007FFFFFU, 0x00FFFFFFU,
    0x01FFFFFFU, 0x03FFFFFFU, 0x07FFFFFFU, 0x0FFFFFFFU,
    0x1FFFFFFFU, 0x3FFFFFFFU, 0x7FFFFFFFU, 0xFFFFFFFFU
};

/* Run length pattern table of 0 (Note: run length from the upper bit) */
const uint32_t g_bitstream_zerobit_runlength_table[256] = {
    8, 7, 6, 6, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4,
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};

/* Table for NLZ calculations */
#define UNUSED 99
static const uint32_t st_nlz10_table[64] = {
        32,     20,     19, UNUSED, UNUSED,     18, UNUSED,      7,
        10,     17, UNUSED, UNUSED,     14, UNUSED,      6, UNUSED,
    UNUSED,      9, UNUSED,     16, UNUSED, UNUSED,      1,     26,
    UNUSED,     13, UNUSED, UNUSED,     24,      5, UNUSED, UNUSED,
    UNUSED,     21, UNUSED,      8,     11, UNUSED,     15, UNUSED,
    UNUSED, UNUSED, UNUSED,      2,     27,      0,     25, UNUSED,
        22, UNUSED,     12, UNUSED, UNUSED,      3,     28, UNUSED,
        23, UNUSED,      4,     29, UNUSED, UNUSED,     30,     31
};
#undef UNUSED

#if !defined(BITSTREAM_USE_MACROS)

/* Open bit reader */
void BitReader_Open(struct BitStream *stream, const uint8_t *memory, size_t size)
{
    /* Argument check */
    assert(stream != NULL);
    assert(memory != NULL);

    /* Reset internal state */
    stream->flags = 0;

    /* Initialize buffer */
    stream->bit_count   = 0;
    stream->bit_buffer  = 0;

    /* memory set */
    stream->memory_image = memory;
    stream->memory_size = size;
    stream->memory_tail = memory + size;

    /* Read position is at the beginning */
    stream->memory_p = (uint8_t *)(memory);

    /* Set as reading mode */
    stream->flags |= (uint8_t)BITSTREAM_FLAGS_MODE_READ;
}

/* Open the bit writer */
void BitWriter_Open(struct BitStream *stream, const uint8_t *memory, size_t size)
{
    /* Argument check */
    assert(stream != NULL);
    assert(memory != NULL);

    /* Reset internal state */
    stream->flags = 0;

    /* Initialize buffer */
    stream->bit_count = 32;
    stream->bit_buffer = 0;

    /* memory set */
    stream->memory_image = memory;
    stream->memory_size = size;
    stream->memory_tail = memory + size;

    /* Read position is at the beginning */
    stream->memory_p = (uint8_t *)(memory);

    /* Set as write mode */
    stream->flags &= (uint8_t)(~BITSTREAM_FLAGS_MODE_READ);
}

/* Close bitstream */
void BitStream_Close(struct BitStream *stream)
{
    /* Argument check */
    assert(stream != NULL);

    /* Flush remaining data */
    BitStream_Flush(stream);

    /* Clear the buffer */
    stream->bit_buffer = 0;

    /* Clear memory information */
    stream->memory_image = NULL;
    stream->memory_size  = 0;

    /* Clear internal state */
    stream->memory_p     = NULL;
    stream->flags        = 0;
}

/* Seek (compliant with fseek) */
void BitStream_Seek(struct BitStream *stream, int32_t offset, int32_t origin)
{
    uint8_t *pos = NULL;

    /* Argument check */
    assert(stream != NULL);

    /* Clear the internal buffer (side effects occur) */
    BitStream_Flush(stream);

    /* First determine the starting point */
    switch (origin) {
    case BITSTREAM_SEEK_CUR:
        pos = stream->memory_p;
        break;
    case BITSTREAM_SEEK_SET:
        pos = (uint8_t *)stream->memory_image;
        break;
    case BITSTREAM_SEEK_END:
        pos = (uint8_t *)((stream)->memory_tail - 1);
        break;
    default:
        assert(0);
    }

    /* Move by offset */
    pos += (offset);

    /* Range check */
    assert(pos >= stream->memory_image);
    assert(pos < (stream)->memory_tail);

    /* Save the results */
    stream->memory_p = pos;
}

/* Based on current position (ftell) */
void BitStream_Tell(struct BitStream *stream, int32_t *result)
{
    /* Argument check */
    assert(stream != NULL);
    assert(result != NULL);

    /* Return the access offset */
    (*result) = (int32_t)(stream->memory_p - stream->memory_image);
}

/* Output n bits to the right (lower) of val (maximum 32 bits can be output) */
void BitWriter_PutBits(struct BitStream *stream, uint32_t val, uint32_t nbits)
{
    /* Argument check */
    assert(stream != NULL);

    /* Not executable in read mode */
    assert(!(stream->flags & BITSTREAM_FLAGS_MODE_READ));

    /* Check if the maximum number of bits that can be output is not exceeded */
    assert(nbits <= 32);

    /* 0 bit output does nothing and ends */
    if (!nbits) { return; }

    /* Outputs val from the most significant bit */
    if (nbits >= stream->bit_count) {
        nbits -= stream->bit_count;
        stream->bit_buffer |= BITSTREAM_GETLOWERBITS(val >> nbits, stream->bit_count);

        /* Check if end is not reached */
        assert(stream->memory_p >= stream->memory_image);
        assert((stream->memory_p + 3) < stream->memory_tail);

        /* Write to memory */
        stream->memory_p[0] = ((stream->bit_buffer >> 24) & 0xFF);
        stream->memory_p[1] = ((stream->bit_buffer >> 16) & 0xFF);
        stream->memory_p[2] = ((stream->bit_buffer >>  8) & 0xFF);
        stream->memory_p[3] = ((stream->bit_buffer >>  0) & 0xFF);
        stream->memory_p += 4;

        /* Reset the buffer */
        stream->bit_buffer = 0;
        stream->bit_count = 32;
    }

    /* Processing fractional bits: Set the remaining bit to the upper bit of the buffer */
    assert(nbits <= 32);
    stream->bit_count -= nbits;
    stream->bit_buffer |= BITSTREAM_GETLOWERBITS(val, nbits) << stream->bit_count;
}

/* Output a run of 0s followed by a terminating 1 */
void BitWriter_PutZeroRun(struct BitStream *stream, uint32_t runlength)
{
    uint32_t run = runlength + 1;

    /* Argument check */
    assert(stream != NULL);

    /* Not executable in read mode */
    assert(!(stream->flags & BITSTREAM_FLAGS_MODE_READ));

    /* Output in 31-bit units */
    while (run > 31) {
        BitWriter_PutBits(stream, 0, 31);
        run -= 31;
    }

    /* Output the terminating 1 */
    BitWriter_PutBits(stream, 1, run);
}

/* Get nbits (up to 32 bits), right-justify the value, and output it */
void BitReader_GetBits(struct BitStream *stream, uint32_t *val, uint32_t nbits)
{
    uint32_t tmp = 0;

    /* Argument check */
    assert(stream != NULL);
    assert(val != NULL);

    /* Assert if not in read mode */
    assert(stream->flags & BITSTREAM_FLAGS_MODE_READ);

    /* Check if the maximum number of bits that can be input is not exceeded */
    assert(nbits <= 32);

    /* Get from buffer */
    if (nbits <= stream->bit_count) {
        stream->bit_count -= nbits;
        (*val) = BITSTREAM_GETLOWERBITS(stream->bit_buffer >> stream->bit_count, nbits);
        return;
    }

    /* If more bits are requested than the current buffer capacity requires, read from memory */

    /* Set the remaining bits to the high order bits */
    nbits -= stream->bit_count;
    tmp = BITSTREAM_GETLOWERBITS(stream->bit_buffer, stream->bit_count) << nbits;

    /* Check if end is not reached */
    assert(stream->memory_p >= stream->memory_image);
    assert(stream->memory_p < stream->memory_tail);

    /* Read from memory */
    stream->bit_buffer
        = ((uint32_t)stream->memory_p[0] << 24) | ((uint32_t)stream->memory_p[1] << 16)
        | ((uint32_t)stream->memory_p[2] << 8) | ((uint32_t)stream->memory_p[3] << 0);
    stream->memory_p += 4;
    stream->bit_count = 32;

    /* Processing of fractional bits Set the remaining bits to the most significant bit of tmp */
    stream->bit_count -= nbits;
    tmp |= BITSTREAM_GETLOWERBITS(stream->bit_buffer >> stream->bit_count, nbits);

    /* normal termination */
    (*val) = tmp;
}

/* Read until the next 1 is encountered, then get the run length of the 0s read in between */
void BitReader_GetZeroRunLength(struct BitStream *stream, uint32_t *runlength)
{
    uint32_t run;

    /* Argument check */
    assert(stream != NULL);
    assert(runlength != NULL);

    /* Measure consecutive 0's from the upper bit */
    run = BITSTREAM_NLZ(BITSTREAM_GETLOWERBITS(stream->bit_buffer, stream->bit_count)) + stream->bit_count - 32;

    /* Decrease the count by the amount read */
    assert(stream->bit_count >= run);
    stream->bit_count -= run;

    /* When the buffer is empty */
    while (!stream->bit_count) {
        /* Read 1 byte and measure again */
        uint32_t tmp_run;

        /* Check if end is not reached */
        assert(stream->memory_p >= stream->memory_image);
        assert(stream->memory_p < stream->memory_tail);

        /* Read from memory, reset to bit buffer and measure run again */
        stream->bit_buffer = stream->memory_p[0];
        stream->memory_p++;
        /* Get run length from table */
        tmp_run = g_bitstream_zerobit_runlength_table[stream->bit_buffer];
        stream->bit_count = 8 - tmp_run;
        /* Add run */
        run += tmp_run;
    }

    /* Read the next 1 blank */
    assert(stream->bit_count >= 1);
    stream->bit_count -= 1;

    /* normal termination */
    (*runlength) = run;
}

/* Clear the bits accumulated in the buffer (move the read/write position to the next byte boundary) */
void BitStream_Flush(struct BitStream *stream)
{
    /* Argument check */
    assert(stream != NULL);

    if (stream->flags & BITSTREAM_FLAGS_MODE_READ) {
        /* Move the read position back by the number of bytes remaining in the buffer and clear the buffer */
        stream->memory_p -= (stream->bit_count >> 3);
        stream->bit_buffer = 0;
        stream->bit_count = 0;
    } else {
        if (stream->bit_count < 32) {
            /* Output up to the next byte boundary */
            const uint32_t remainbits = 32 - stream->bit_count;
            if (remainbits > 24) {
                stream->memory_p[0] = ((stream->bit_buffer >> 24) & 0xFF);
                stream->memory_p[1] = ((stream->bit_buffer >> 16) & 0xFF);
                stream->memory_p[2] = ((stream->bit_buffer >>  8) & 0xFF);
                stream->memory_p[3] = ((stream->bit_buffer >>  0) & 0xFF);
                stream->memory_p += 4;
            } else if (remainbits > 16) {
                stream->memory_p[0] = ((stream->bit_buffer >> 24) & 0xFF);
                stream->memory_p[1] = ((stream->bit_buffer >> 16) & 0xFF);
                stream->memory_p[2] = ((stream->bit_buffer >>  8) & 0xFF);
                stream->memory_p += 3;
            } else if (remainbits > 8) {
                stream->memory_p[0] = ((stream->bit_buffer >> 24) & 0xFF);
                stream->memory_p[1] = ((stream->bit_buffer >> 16) & 0xFF);
                stream->memory_p += 2;
            } else {
                stream->memory_p[0] = ((stream->bit_buffer >> 24) & 0xFF);
                stream->memory_p += 1;
            }
            stream->bit_count = 32;
            stream->bit_buffer = 0;
        }
    }
}

#endif /* BITSTREAM_USE_MACROS */

/* Calculate NLZ (number of bits from the most significant bit to 1) */
uint32_t BitStream_NLZSoft(uint32_t x)
{
    /* Hacker's Fun Reference */
    x = x | (x >> 1);
    x = x | (x >> 2);
    x = x | (x >> 4);
    x = x | (x >> 8);
    x = x & ~(x >> 16);
    x = (x << 9) - x;
    x = (x << 11) - x;
    x = (x << 14) - x;
    return st_nlz10_table[x >> 26];
}
