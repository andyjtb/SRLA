#ifndef BYTEARRAY_H_INCLUDED
#define BYTEARRAY_H_INCLUDED

#include <stdint.h>

/* Read 1 byte */
#define ByteArray_ReadUint8(p_array)\
    (uint8_t)((p_array)[0])

/* Read 2 bytes (big endian) */
#define ByteArray_ReadUint16BE(p_array)\
    (uint16_t)(\
            (((uint16_t)((p_array)[0])) << 8) |\
            (((uint16_t)((p_array)[1])) << 0)\
            )

/* Read 3 bytes (big endian) */
#define ByteArray_ReadUint24BE(p_array)\
    (uint32_t)(\
            (((uint32_t)((p_array)[0])) << 16) |\
            (((uint32_t)((p_array)[1])) <<  8) |\
            (((uint32_t)((p_array)[2])) <<  0)\
            )

/* Read 4 bytes (big endian) */
#define ByteArray_ReadUint32BE(p_array)\
    (uint32_t)(\
            (((uint32_t)((p_array)[0])) << 24) |\
            (((uint32_t)((p_array)[1])) << 16) |\
            (((uint32_t)((p_array)[2])) <<  8) |\
            (((uint32_t)((p_array)[3])) <<  0)\
            )

/* Read 2 bytes (little endian) */
#define ByteArray_ReadUint16LE(p_array)\
    (uint16_t)(\
            (((uint16_t)((p_array)[0])) << 0) |\
            (((uint16_t)((p_array)[1])) << 8)\
            )

/* Read 3 bytes (little endian) */
#define ByteArray_ReadUint24LE(p_array)\
    (uint32_t)(\
            (((uint32_t)((p_array)[0])) <<  0) |\
            (((uint32_t)((p_array)[1])) <<  8) |\
            (((uint32_t)((p_array)[2])) << 16)\
            )

/* Read 4 bytes (little endian) */
#define ByteArray_ReadUint32LE(p_array)\
    (uint32_t)(\
            (((uint32_t)((p_array)[0])) <<  0) |\
            (((uint32_t)((p_array)[1])) <<  8) |\
            (((uint32_t)((p_array)[2])) << 16) |\
            (((uint32_t)((p_array)[3])) << 24)\
            )

/* Get 1 byte */
#define ByteArray_GetUint8(p_array, p_u8val)\
    do {\
        (*(p_u8val)) = ByteArray_ReadUint8(p_array);\
        (p_array) += 1;\
    } while (0);

/* Get 2 bytes (big endian) */
#define ByteArray_GetUint16BE(p_array, p_u16val)\
    do {\
        (*(p_u16val)) = ByteArray_ReadUint16BE(p_array);\
        (p_array) += 2;\
    } while (0);

/* Get 3 bytes (big endian) */
#define ByteArray_GetUint24BE(p_array, p_u32val)\
    do {\
        (*(p_u32val)) = ByteArray_ReadUint24BE(p_array);\
        (p_array) += 3;\
    } while (0);

/* Get 4 bytes (big endian) */
#define ByteArray_GetUint32BE(p_array, p_u32val)\
    do {\
        (*(p_u32val)) = ByteArray_ReadUint32BE(p_array);\
        (p_array) += 4;\
    } while (0);

/* Get 2 bytes (little endian) */
#define ByteArray_GetUint16LE(p_array, p_u16val)\
    do {\
        (*(p_u16val)) = ByteArray_ReadUint16LE(p_array);\
        (p_array) += 2;\
    } while (0);

/* Get 3 bytes (little endian) */
#define ByteArray_GetUint24LE(p_array, p_u32val)\
    do {\
        (*(p_u32val)) = ByteArray_ReadUint24LE(p_array);\
        (p_array) += 3;\
    } while (0);

/* Get 4 bytes (little endian) */
#define ByteArray_GetUint32LE(p_array, p_u32val)\
    do {\
        (*(p_u32val)) = ByteArray_ReadUint32LE(p_array);\
        (p_array) += 4;\
    } while (0);

/* Write 1 byte */
#define ByteArray_WriteUint8(p_array, u8val)\
    do {\
        ((p_array)[0]) = (uint8_t)(u8val);\
    } while (0);

/* Write 2 bytes (big endian) */
#define ByteArray_WriteUint16BE(p_array, u16val)\
    do {\
        ((p_array)[0]) = (uint8_t)(((u16val) >> 8) & 0xFF);\
        ((p_array)[1]) = (uint8_t)(((u16val) >> 0) & 0xFF);\
    } while (0);

/* Write 3 bytes (big endian) */
#define ByteArray_WriteUint24BE(p_array, u32val)\
    do {\
        ((p_array)[0]) = (uint8_t)(((u32val) >> 16) & 0xFF);\
        ((p_array)[1]) = (uint8_t)(((u32val) >>  8) & 0xFF);\
        ((p_array)[2]) = (uint8_t)(((u32val) >>  0) & 0xFF);\
    } while (0);

/* Write 4 bytes (big endian) */
#define ByteArray_WriteUint32BE(p_array, u32val)\
    do {\
        ((p_array)[0]) = (uint8_t)(((u32val) >> 24) & 0xFF);\
        ((p_array)[1]) = (uint8_t)(((u32val) >> 16) & 0xFF);\
        ((p_array)[2]) = (uint8_t)(((u32val) >>  8) & 0xFF);\
        ((p_array)[3]) = (uint8_t)(((u32val) >>  0) & 0xFF);\
    } while (0);

/* Write 2 bytes (little endian) */
#define ByteArray_WriteUint16LE(p_array, u16val)\
    do {\
        ((p_array)[0]) = (uint8_t)(((u16val) >> 0) & 0xFF);\
        ((p_array)[1]) = (uint8_t)(((u16val) >> 8) & 0xFF);\
    } while (0);

/* Write 3 bytes (little endian) */
#define ByteArray_WriteUint24LE(p_array, u32val)\
    do {\
        ((p_array)[0]) = (uint8_t)(((u32val) >>  0) & 0xFF);\
        ((p_array)[1]) = (uint8_t)(((u32val) >>  8) & 0xFF);\
        ((p_array)[2]) = (uint8_t)(((u32val) >> 16) & 0xFF);\
    } while (0);

/* Write 4 bytes (little endian) */
#define ByteArray_WriteUint32LE(p_array, u32val)\
    do {\
        ((p_array)[0]) = (uint8_t)(((u32val) >>  0) & 0xFF);\
        ((p_array)[1]) = (uint8_t)(((u32val) >>  8) & 0xFF);\
        ((p_array)[2]) = (uint8_t)(((u32val) >> 16) & 0xFF);\
        ((p_array)[3]) = (uint8_t)(((u32val) >> 24) & 0xFF);\
    } while (0);

/* 1 byte output */
#define ByteArray_PutUint8(p_array, u8val)\
    do {\
        ByteArray_WriteUint8(p_array, u8val);\
        (p_array) += 1;\
    } while (0);

/* 2-byte output (big endian) */
#define ByteArray_PutUint16BE(p_array, u16val)\
    do {\
        ByteArray_WriteUint16BE(p_array, u16val);\
        (p_array) += 2;\
    } while (0);

/* 3-byte output (big endian) */
#define ByteArray_PutUint24BE(p_array, u32val)\
    do {\
        ByteArray_WriteUint24BE(p_array, u32val);\
        (p_array) += 3;\
    } while (0);

/* 4-byte output (big endian) */
#define ByteArray_PutUint32BE(p_array, u32val)\
    do {\
        ByteArray_WriteUint32BE(p_array, u32val);\
        (p_array) += 4;\
    } while (0);

/* 2-byte output (little endian) */
#define ByteArray_PutUint16LE(p_array, u16val)\
    do {\
        ByteArray_WriteUint16LE(p_array, u16val);\
        (p_array) += 2;\
    } while (0);

/* 3-byte output (little endian) */
#define ByteArray_PutUint24LE(p_array, u32val)\
    do {\
        ByteArray_WriteUint24LE(p_array, u32val);\
        (p_array) += 3;\
    } while (0);

/* 4-byte output (little endian) */
#define ByteArray_PutUint32LE(p_array, u32val)\
    do {\
        ByteArray_WriteUint32LE(p_array, u32val);\
        (p_array) += 4;\
    } while (0);

#endif /* BYTEARRAY_H_INCLUDED */
