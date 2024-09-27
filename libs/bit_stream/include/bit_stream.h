#ifndef BITSTREAM_H_INCLUDED
#define BITSTREAM_H_INCLUDED

#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <limits.h>

#if CHAR_BIT != 8
#error "This program 8bit/byte system only."
#endif

/* Use macros or not? */
#define BITSTREAM_USE_MACROS 1

/* Search code for BitStream_Seek function */
#define BITSTREAM_SEEK_SET  (int32_t)SEEK_SET
#define BITSTREAM_SEEK_CUR  (int32_t)SEEK_CUR
#define BITSTREAM_SEEK_END  (int32_t)SEEK_END

/* Read mode? (0 is write mode) */
#define BITSTREAM_FLAGS_MODE_READ  (1 << 0)

/* Get the lowest nbits of val */
#define BITSTREAM_GETLOWERBITS(val, nbits) ((val) & g_bitstream_lower_bits_mask[(nbits)])

/* Bitstream structure */
struct BitStream {
    uint32_t bit_buffer; /* temporary buffer of bits */
    uint32_t bit_count; /* [Reader] Number of remaining buffer bits, [Writer] Number of bits until memory write */
    const uint8_t *memory_image; /* Start of memory area */
    const uint8_t *memory_tail;/* End of memory area */
    size_t memory_size; /* Memory area size */
    uint8_t *memory_p; /* Memory read/write location */
    uint8_t flags; /* Internal state flag */
};

/* Calculate NLZ (number of bits from the most significant bit to 1) */
#if defined(__GNUC__)
/* Use built-in functions */
#define BITSTREAM_NLZ(x) (((x) > 0) ? (uint32_t)__builtin_clz(x) : 32U)
#elif defined(_MSC_VER)
#include <intrin.h>
/* Use built-in functions */
__inline uint32_t BITSTREAM_NLZ(uint32_t x)
{
    return __lzcnt(x);
}
#else
/* Use software implementation */
#define BITSTREAM_NLZ(x) BitStream_NLZSoft(x)
#endif

/* Software implementation of NLZ */
#ifdef __cplusplus
extern "C" uint32_t BitStream_NLZSoft(uint32_t x);
#else
uint32_t BitStream_NLZSoft(uint32_t x);
#endif

#if !defined(BITSTREAM_USE_MACROS)

#ifdef __cplusplus
extern "C" {
#endif

/* Open bit reader */
void BitReader_Open(struct BitStream *stream, const uint8_t *memory, size_t size);

/* Open the bit writer */
void BitWriter_Open(struct BitStream *stream, const uint8_t *memory, size_t size);

/* Close bitstream */
void BitStream_Close(struct BitStream *stream);

/* Seek (compliant with fseek) */
void BitStream_Seek(struct BitStream *stream, int32_t offset, int32_t origin);

/* Based on current position (ftell) */
void BitStream_Tell(struct BitStream *stream, int32_t *result);

/* Output n bits to the right (lower) of val (maximum 32 bits can be output) */
void BitWriter_PutBits(struct BitStream *stream, uint32_t val, uint32_t nbits);

/* Output a run of 0s followed by a terminating 1 */
void BitWriter_PutZeroRun(struct BitStream *stream, uint32_t runlength);

/* Get nbits (up to 32 bits), right-justify the value, and output it */
void BitReader_GetBits(struct BitStream *stream, uint32_t *val, uint32_t nbits);

/* Read until the next 1 is encountered, then get the run length of the 0s read in between */
void BitReader_GetZeroRunLength(struct BitStream *stream, uint32_t *runlength);

/* Clear the bits accumulated in the buffer */
void BitStream_Flush(struct BitStream *stream);

#ifdef __cplusplus
}
#endif

#else /* BITSTREAM_USE_MACROS */

#ifdef __cplusplus
extern "C" {
#endif

/* Mask to extract the lower bits */
extern const uint32_t g_bitstream_lower_bits_mask[33];

/* Run length pattern table */
extern const uint32_t g_bitstream_zerobit_runlength_table[0x100];

#ifdef __cplusplus
}
#endif

/* Open bit reader */
#define BitReader_Open(stream, memory, size)\
    do {\
        /* Argument check */\
        assert((void *)(stream) != NULL);\
        assert((void *)(memory) != NULL);\
        \
        /* Reset internal state */\
        (stream)->flags = 0;\
        \
        /* Initialize buffer */\
        (stream)->bit_count   = 0;\
        (stream)->bit_buffer  = 0;\
        \
        /* memory set */\
        (stream)->memory_image = (memory);\
        (stream)->memory_size  = (size);\
        (stream)->memory_tail  = (memory) + (size);\
        \
        /* Read position is at the beginning */\
        (stream)->memory_p = (memory);\
        \
        /* Set as reading mode */\
        (stream)->flags |= (uint8_t)BITSTREAM_FLAGS_MODE_READ;\
    } while (0)

/* Open the bit writer */
#define BitWriter_Open(stream, memory, size)\
    do {\
        /* Argument check */\
        assert((void *)(stream) != NULL);\
        assert((void *)(memory) != NULL);\
        \
        /* Reset internal state */\
        (stream)->flags = 0;\
        \
        /* Initialize buffer */\
        (stream)->bit_count   = 32;\
        (stream)->bit_buffer  = 0;\
        \
        /* memory set */\
        (stream)->memory_image = (memory);\
        (stream)->memory_size  = (size);\
        (stream)->memory_tail  = (memory) + (size);\
        \
        /* Read position is at the beginning */\
        (stream)->memory_p = (memory);\
        \
        /* Set as write mode */\
        (stream)->flags &= (uint8_t)(~BITSTREAM_FLAGS_MODE_READ);\
    } while (0)

/* Close bitstream */
#define BitStream_Close(stream)\
    do {\
        /* Argument check */\
        assert((void *)(stream) != NULL);\
        \
        /* Flush remaining data */\
        BitStream_Flush(stream);\
        \
        /* Clear the buffer */\
        (stream)->bit_buffer = 0;\
        \
        /* Clear memory information */\
        (stream)->memory_image = NULL;\
        (stream)->memory_size  = 0;\
        \
        /* Clear internal state */\
        (stream)->memory_p     = NULL;\
        (stream)->flags        = 0;\
    } while (0)

/* Seek (compliant with fseek) */
#define BitStream_Seek(stream, offset, origin)\
    do {\
        uint8_t* __pos = NULL;\
        \
        /* Argument check */\
        assert((void *)(stream) != NULL);\
        \
        /* Clear the internal buffer (side effects occur) */\
        BitStream_Flush(stream);\
        \
        /* First determine the starting point */\
        switch (origin) {\
        case BITSTREAM_SEEK_CUR:\
            __pos = (stream)->memory_p;\
            break;\
        case BITSTREAM_SEEK_SET:\
            __pos = (uint8_t *)(stream)->memory_image;\
            break;\
        case BITSTREAM_SEEK_END:\
            __pos = (uint8_t *)((stream)->memory_tail - 1);\
            break;\
        default:\
            assert(0);\
        }\
        \
        /* Move by offset */\
        __pos += (offset);\
        \
        /* Range check */\
        assert(__pos >= (stream)->memory_image);\
        assert(__pos < (stream)->memory_tail);\
        \
        /* Save the results */\
        (stream)->memory_p = __pos;\
    } while (0)

/* Based on current position (ftell) */
#define BitStream_Tell(stream, result)\
    do {\
        /* Argument check */\
        assert((void *)(stream) != NULL);\
        assert((void *)(result) != NULL);\
        \
        /* Return the access offset */\
        (*result) = (int32_t)\
        ((stream)->memory_p - (stream)->memory_image);\
    } while (0)

/* Output n bits to the right (lower) of val (maximum 32 bits can be output) */
#define BitWriter_PutBits(stream, val, nbits)\
    do {\
        uint32_t __nbits;\
        \
        /* Argument check */\
        assert((void *)(stream) != NULL);\
        \
        /* Not executable in read mode */\
        assert(!((stream)->flags & BITSTREAM_FLAGS_MODE_READ));\
        \
        /* Check if the maximum number of bits that can be output is not exceeded */\
        assert((nbits) <= 32);\
        \
        /* 0 bit output does nothing and ends */\
        if (!(nbits)) { break; }\
        \
        /* Outputs val from the most significant bit */\
        __nbits = (nbits);\
        if (__nbits >= (stream)->bit_count) {\
            __nbits -= (stream)->bit_count;\
            (stream)->bit_buffer |= BITSTREAM_GETLOWERBITS((val) >> __nbits, (stream)->bit_count);\
            \
            /* Check if end is not reached */\
            assert((stream)->memory_p >= (stream)->memory_image);\
            assert(((stream)->memory_p + 3) < (stream)->memory_tail);\
            \
            /* Write to memory */\
            (stream)->memory_p[0] = (uint8_t)(((stream)->bit_buffer >> 24) & 0xFF);\
            (stream)->memory_p[1] = (uint8_t)(((stream)->bit_buffer >> 16) & 0xFF);\
            (stream)->memory_p[2] = (uint8_t)(((stream)->bit_buffer >>  8) & 0xFF);\
            (stream)->memory_p[3] = (uint8_t)(((stream)->bit_buffer >>  0) & 0xFF);\
            (stream)->memory_p += 4;\
            \
            /* Reset the buffer */\
            (stream)->bit_buffer = 0;\
            (stream)->bit_count = 32;\
        }\
        \
        /* Processing fractional bits: Set the remaining bit to the upper bit of the buffer */\
        assert(__nbits <= 32);\
        (stream)->bit_count -= __nbits;\
        (stream)->bit_buffer |= BITSTREAM_GETLOWERBITS(val, __nbits) << (stream)->bit_count;\
    } while (0)

/* Output a run of 0s followed by a terminating 1 */
#define BitWriter_PutZeroRun(stream, runlength)\
    do {\
        uint32_t __run = ((runlength) + 1);\
        \
        /* Argument check */\
        assert((void *)(stream) != NULL);\
        \
        /* Not executable in read mode */\
        assert(!((stream)->flags & BITSTREAM_FLAGS_MODE_READ));\
        \
        /* Output in 31-bit units */\
        while (__run > 31) {\
            BitWriter_PutBits(stream, 0, 31);\
            __run -= 31;\
        }\
        /* Output the terminating 1 */\
        BitWriter_PutBits(stream, 1, __run);\
    } while (0)

/* Get nbits (up to 32 bits), right-justify the value, and output it */
#define BitReader_GetBits(stream, val, nbits)\
    do {\
        uint32_t __tmp, __nbits;\
        \
        /* Argument check */\
        assert((void *)(stream) != NULL);\
        assert((void *)(val) != NULL);\
        \
        /* Assert if not in read mode */\
        assert((stream)->flags & BITSTREAM_FLAGS_MODE_READ);\
        \
        /* Check if the maximum number of bits that can be input is not exceeded */\
        assert((nbits) <= 32);\
        \
        /* Get from buffer */\
        if ((nbits) <= (stream)->bit_count) {\
            (stream)->bit_count -= (nbits);\
            (*(val)) = BITSTREAM_GETLOWERBITS((stream)->bit_buffer >> (stream)->bit_count, (nbits));\
            break;\
        }\
        \
        /* If more bits are requested than the current buffer capacity requires, read from memory */\
        __nbits = (nbits);\
        /* Set the remaining bits to the high order bits */\
        __nbits -= (stream)->bit_count;\
        __tmp = BITSTREAM_GETLOWERBITS((stream)->bit_buffer, (stream)->bit_count) << __nbits;\
        \
        /* Check if end is not reached */\
        assert((stream)->memory_p >= (stream)->memory_image);\
        assert((stream)->memory_p < (stream)->memory_tail);\
        \
        /* Read from memory */\
        (stream)->bit_buffer\
            = ((uint32_t)(stream)->memory_p[0] << 24) | ((uint32_t)(stream)->memory_p[1] << 16)\
            | ((uint32_t)(stream)->memory_p[2] <<  8) | ((uint32_t)(stream)->memory_p[3] <<  0);\
        (stream)->memory_p += 4;\
        (stream)->bit_count = 32;\
        \
        /* Processing of fractional bits Set the remaining bits to the most significant bit of tmp */\
        (stream)->bit_count -= __nbits;\
        __tmp |= BITSTREAM_GETLOWERBITS((stream)->bit_buffer >> (stream)->bit_count, __nbits);\
        \
        /* normal termination */\
        (*(val)) = __tmp;\
    } while (0)

/* Read until the next 1 is encountered, then get the run length of the 0s read in between */
#define BitReader_GetZeroRunLength(stream, runlength)\
    do {\
        uint32_t __run;\
        \
        /* Argument check */\
        assert((void *)(stream) != NULL);\
        assert((void *)(runlength) != NULL);\
        \
        /* Measure consecutive 0's from the upper bit */\
        __run = BITSTREAM_NLZ(BITSTREAM_GETLOWERBITS((stream)->bit_buffer, (stream)->bit_count)) + (stream)->bit_count - 32;\
        \
        /* Decrease the count by the amount read */\
        assert((stream)->bit_count >= __run);\
        (stream)->bit_count -= __run;\
        \
        /* When the buffer is empty */\
        while (!(stream)->bit_count) {\
            /* Read 1 byte and measure again */\
            uint32_t __tmp_run;\
            \
            /* Check if end is not reached */\
            assert((stream)->memory_p >= (stream)->memory_image);\
            assert((stream)->memory_p < (stream)->memory_tail);\
            \
            /* Read from memory, reset to bit buffer and measure run again */\
            (stream)->bit_buffer = (stream)->memory_p[0];\
            (stream)->memory_p++;\
            /* Get run length from table */\
            __tmp_run = g_bitstream_zerobit_runlength_table[(stream)->bit_buffer];\
            (stream)->bit_count = 8 - __tmp_run;\
            /* Add run */\
            __run += __tmp_run;\
        }\
        \
        /* Read the next 1 blank */\
        assert((stream)->bit_count >= 1);\
        (stream)->bit_count -= 1;\
        \
        /* normal termination */\
        (*(runlength)) = __run;\
    } while (0)

/* Clear the bits accumulated in the buffer */
#define BitStream_Flush(stream)\
    do {\
        /* Argument check */\
        assert((void *)(stream) != NULL);\
        \
        if ((stream)->flags & BITSTREAM_FLAGS_MODE_READ) {\
            /* Move the read position back by the number of bytes remaining in the buffer and clear the buffer */\
            (stream)->memory_p -= ((stream)->bit_count >> 3);\
            (stream)->bit_buffer = 0;\
            (stream)->bit_count = 0;\
        } else {\
            if ((stream)->bit_count < 32) {\
                /* Output up to the next byte boundary */\
                const uint32_t __remainbits = 32 - (stream)->bit_count;\
                if (__remainbits > 24) {\
                    (stream)->memory_p[0] = (uint8_t)(((stream)->bit_buffer >> 24) & 0xFF);\
                    (stream)->memory_p[1] = (uint8_t)(((stream)->bit_buffer >> 16) & 0xFF);\
                    (stream)->memory_p[2] = (uint8_t)(((stream)->bit_buffer >>  8) & 0xFF);\
                    (stream)->memory_p[3] = (uint8_t)(((stream)->bit_buffer >>  0) & 0xFF);\
                    (stream)->memory_p += 4;\
                } else if (__remainbits > 16) {\
                    (stream)->memory_p[0] = (uint8_t)(((stream)->bit_buffer >> 24) & 0xFF);\
                    (stream)->memory_p[1] = (uint8_t)(((stream)->bit_buffer >> 16) & 0xFF);\
                    (stream)->memory_p[2] = (uint8_t)(((stream)->bit_buffer >>  8) & 0xFF);\
                    (stream)->memory_p += 3;\
                } else if (__remainbits > 8) {\
                    (stream)->memory_p[0] = (uint8_t)(((stream)->bit_buffer >> 24) & 0xFF);\
                    (stream)->memory_p[1] = (uint8_t)(((stream)->bit_buffer >> 16) & 0xFF);\
                    (stream)->memory_p += 2;\
                } else {\
                    (stream)->memory_p[0] = (uint8_t)(((stream)->bit_buffer >> 24) & 0xFF);\
                    (stream)->memory_p += 1;\
                }\
                (stream)->bit_count = 32;\
                (stream)->bit_buffer = 0;\
            }\
        }\
    } while (0)

#endif /* BITSTREAM_USE_MACROS */

#endif /* BITSTREAM_H_INCLUDED */
