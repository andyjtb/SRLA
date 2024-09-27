#include "srla_player.h"
#include <assert.h>

#include <AudioToolbox/AudioQueue.h>
#include <CoreAudio/CoreAudioTypes.h>
#include <CoreFoundation/CFRunLoop.h>

#define NUM_BUFFERS 3
#define BUFFER_SIZE (8 * 1024)
#define DECODE_BUFFER_NUM_SAMPLES (1024)

/* CoreAudio output callback function */
static void SRLAPlayer_CoreAudioCallback(void *custom_data, AudioQueueRef queue, AudioQueueBufferRef buffer);

/* Initialize count */
static int32_t st_initialize_count = 0;
/* Player config at initialization */
static struct SRLAPlayerConfig st_config = { 0, };
/* Buffer area for decoded data */
static int32_t **st_decode_buffer = NULL;
/* Buffer reference position */
static uint32_t st_buffer_pos = DECODE_BUFFER_NUM_SAMPLES; /* Empty state */
/* Reference to output queue */
static AudioQueueRef queue = NULL;

/* Initialization This function initializes the device driver and starts playback. */
void SRLAPlayer_Initialize(const struct SRLAPlayerConfig *config)
{
    uint32_t i;
    AudioStreamBasicDescription format;
    AudioQueueBufferRef buffers[NUM_BUFFERS];

    assert(config != NULL);

    /* Multiple initialization is not possible */
    if (st_initialize_count > 0) {
        return;
    }

    /* Get config */
    st_config = (*config);

    /* Fill the format with attributes */
    format.mSampleRate       = st_config.sampling_rate; /* Sampling rate */
    format.mFormatID         = kAudioFormatLinearPCM; /* Format: PCM */
    format.mFormatFlags      = kLinearPCMFormatFlagIsSignedInteger | kAudioFormatFlagIsPacked; /* Specify format flags. */
    format.mBitsPerChannel   = st_config.bits_per_sample; /* Number of bits per channel */
    format.mChannelsPerFrame = st_config.num_channels; /* Number of channels */
    format.mBytesPerFrame    = (st_config.bits_per_sample * st_config.num_channels) / 8; /* Number of bytes in one frame (one sample of all frames) */
    format.mFramesPerPacket  = 1; /* Number of frames per packet */
    format.mBytesPerPacket   = format.mBytesPerFrame * format.mFramesPerPacket; /* Number of bytes per packet */
    format.mReserved         = 0; /* (reserved area) */

    /* Allocate buffer for decoding area */
    st_decode_buffer = (int32_t **)malloc(sizeof(int32_t *) * st_config.num_channels);
    for (i = 0; i < st_config.num_channels; i++) {
        st_decode_buffer[i] = (int32_t *)malloc(sizeof(int32_t) * DECODE_BUFFER_NUM_SAMPLES);
        memset(st_decode_buffer[i], 0, sizeof(int32_t) * DECODE_BUFFER_NUM_SAMPLES);
    }

    /* Create a new output queue */
    AudioQueueNewOutput(&format,
            SRLAPlayer_CoreAudioCallback, NULL, CFRunLoopGetCurrent(), kCFRunLoopCommonModes, 0, &queue);

    for (i = 0; i < NUM_BUFFERS; i++) {
        /* Allocate buffer space for the specified queue */
        AudioQueueAllocateBuffer(queue, BUFFER_SIZE, &buffers[i]);
        /* Set the size */
        buffers[i]->mAudioDataByteSize = BUFFER_SIZE;
        /* Output the first data */
        SRLAPlayer_CoreAudioCallback(NULL, queue, buffers[i]);
    }

    /* Start playing the cue */
    AudioQueueStart(queue, NULL);

    /*
/* Start of thread loop processing
* The thread processing will continue even after this function ends (apparently called a monitor loop) */
*/
    CFRunLoopRun();

    st_initialize_count++;
}

/* End. Release the initialized resources here. */
void SRLAPlayer_Finalize(void)
{
    if (st_initialize_count == 1) {
        uint32_t i;
        /* Stop and discard the queue */
        AudioQueueStop(queue, false);
        AudioQueueDispose(queue, false);
        CFRunLoopStop(CFRunLoopGetCurrent());

        /* Free the buffer for the decoding area */
        for (i = 0; i < st_config.num_channels; i++) {
            free(st_decode_buffer[i]);
        }
        free(st_decode_buffer);
    }

    st_initialize_count--;
}

/* CoreAudio output callback function */
static void SRLAPlayer_CoreAudioCallback(void *custom_data, AudioQueueRef queue, AudioQueueBufferRef buffer)
{
    uint32_t i, ch;
    int16_t *ch_interleaved_buffer = (int16_t *)buffer->mAudioData;
    const uint32_t num_buffer_samples = BUFFER_SIZE / sizeof(int16_t);

    for (i = 0; i < num_buffer_samples; i += st_config.num_channels) {
        /* If the buffer is full, request the next data immediately */
        if (st_buffer_pos >= DECODE_BUFFER_NUM_SAMPLES) {
            st_config.sample_request_callback(st_decode_buffer, st_config.num_channels, DECODE_BUFFER_NUM_SAMPLES);
            st_buffer_pos = 0;
        }

        /* Fill the interleaved buffer with data */
        for (ch = 0; ch < st_config.num_channels; ch++) {
            ch_interleaved_buffer[i + ch] = (int16_t)st_decode_buffer[ch][st_buffer_pos];
        }
        st_buffer_pos++;
    }

    /* Enqueue the buffer */
    AudioQueueEnqueueBuffer(queue, buffer, 0, NULL);
}
