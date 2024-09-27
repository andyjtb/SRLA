#include "srla_player.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

/* Disable the inline keyword to build in a C89 environment */
#define inline
#include <pulse/simple.h>
#include <pulse/error.h>
#include <pulse/gccmacro.h>

#define BUFFER_SIZE (1 * 1024)
#define DECODE_BUFFER_NUM_SAMPLES (1024)

/* Initialize count */
static int32_t st_initialize_count = 0;
/* Player config at initialization */
static struct SRLAPlayerConfig st_config = { 0, };
/* Buffer area for decoded data */
static int32_t **st_decode_buffer = NULL;
/* Buffer reference position */
static uint32_t st_buffer_pos = DECODE_BUFFER_NUM_SAMPLES; /* Empty state */
/* simple pulseaudio handle */
static pa_simple *pa_simple_hn = NULL;
/* Buffer area */
static uint8_t buffer[BUFFER_SIZE];

/* Initialization This function initializes the device driver and starts playback. */
void SRLAPlayer_Initialize(const struct SRLAPlayerConfig *config)
{
    uint32_t i;
    int error;
    pa_sample_spec sample_spec;

    assert(config != NULL);

    /* Multiple initialization is not possible */
    if (st_initialize_count > 0) {
        return;
    }

    /* Get config */
    st_config = (*config);

    /* Fill the format with attributes */
    sample_spec.format = PA_SAMPLE_S16LE;
    sample_spec.rate = st_config.sampling_rate;
    sample_spec.channels = (uint8_t)st_config.num_channels;

    /* Allocate buffer for decoding area */
    st_decode_buffer = (int32_t **)malloc(sizeof(int32_t *) * st_config.num_channels);
    for (i = 0; i < st_config.num_channels; i++) {
        st_decode_buffer[i] = (int32_t *)malloc(sizeof(int32_t) * DECODE_BUFFER_NUM_SAMPLES);
        memset(st_decode_buffer[i], 0, sizeof(int32_t) * DECODE_BUFFER_NUM_SAMPLES);
    }

    /* Create playback handle */
    if ((pa_simple_hn = pa_simple_new(NULL, "SRLAPlayer", PA_STREAM_PLAYBACK, NULL, "playback",
                    &sample_spec, NULL, NULL, &error)) == NULL) {
        fprintf(stderr, "failed to create pulseaudio playback: %s \n", pa_strerror(error));
        exit(1);
    }

    st_initialize_count++;

    while (1) {
        uint32_t i, ch;
        int16_t *pbuffer = (int16_t *)&buffer[0];
        const uint32_t num_writable_samples_per_channel = (uint32_t)(BUFFER_SIZE / (st_config.num_channels * sizeof(int16_t)));

        for (i = 0; i < num_writable_samples_per_channel; i++) {
            /* If the buffer is full, request the next data immediately */
            if (st_buffer_pos >= DECODE_BUFFER_NUM_SAMPLES) {
                st_config.sample_request_callback(st_decode_buffer, st_config.num_channels, DECODE_BUFFER_NUM_SAMPLES);
                st_buffer_pos = 0;
            }
            /* Fill the interleaved buffer with data */
            for (ch = 0; ch < st_config.num_channels; ch++) {
                *pbuffer++ = (int16_t)st_decode_buffer[ch][st_buffer_pos];
            }
            st_buffer_pos++;
        }

        if (pa_simple_write(pa_simple_hn, buffer, BUFFER_SIZE, &error) < 0) {
            fprintf(stderr, "pa_simple_write() failed: %s\n", pa_strerror(error));
            exit(1);
        }
    }
}

/* End. Release the initialized resources here. */
void SRLAPlayer_Finalize(void)
{
    if (st_initialize_count == 1) {
        uint32_t i;

        pa_simple_free(pa_simple_hn);

        /* Free the buffer for the decoding area */
        for (i = 0; i < st_config.num_channels; i++) {
            free(st_decode_buffer[i]);
        }
        free(st_decode_buffer);
    }

    st_initialize_count--;
}
