#include "srla_player.h"
#include <assert.h>
#include <stdio.h>

#include <windows.h>

/* Define COBJMACROS to use macro calls */
#define COBJMACROS
#include <mmdeviceapi.h>
#include <audioclient.h>
#undef COBJMACROS

#define DECODE_BUFFER_NUM_SAMPLES (1024)
#define REQUESTED_SOUND_BUFFER_DURATION  (2 * 10000000LL) /* Internally requested buffer size [100 nanoseconds] */

/* Initialize count */
static int32_t st_initialize_count = 0;
/* Player config at initialization */
static struct SRLAPlayerConfig st_config = { 0, };
/* Buffer area for decoded data */
static int32_t** st_decode_buffer = NULL;
/* Buffer reference position */
static uint32_t st_buffer_pos = DECODE_BUFFER_NUM_SAMPLES; /* Empty state */
/* Handle for WASAPI control */
static IAudioClient* audio_client = NULL;
static IAudioRenderClient* audio_render_client = NULL;

/* Define your own CLSID and IID */
/* Supplementary Note) __uuid cannot be used unless the source is C++. If you want to use C++, use classes. However, I don't want to make everything a C++ project just for Windows reasons. */
static const CLSID st_CLSID_MMDeviceEnumerator = { 0xBCDE0395, 0xE52F, 0x467C, {0x8E,0x3D,0xC4,0x57,0x92,0x91,0x69,0x2E} };
static const IID st_IID_IMMDeviceEnumerator = { 0xA95664D2, 0x9614, 0x4F35, {0xA7,0x46,0xDE,0x8D,0xB6,0x36,0x17,0xE6} };
static const IID st_IID_IAudioClient = { 0x1CB9AD4C, 0xDBFA, 0x4C32, {0xB1,0x78,0xC2,0xF5,0x68,0xA7,0x03,0xB2} };
static const IID st_IID_IAudioClockAdjustment = { 0xF6E4C0A0, 0x46D9, 0x4FB8, {0xBE,0x21,0x57,0xA3,0xEF,0x2B,0x62,0x6C} };
static const IID st_IID_IAudioRenderClient = { 0xF294ACFC, 0x3146, 0x4483, {0xA7,0xBF,0xAD,0xDC,0xA7,0xC2,0x60,0xE2} };

/* Initialization This function initializes the device driver and starts playback. */
void SRLAPlayer_Initialize(const struct SRLAPlayerConfig* config)
{
    uint32_t i, buffer_frame_size;
    HRESULT hr;
    IMMDeviceEnumerator* device_enumerator;
    IMMDevice* audio_device;
    WAVEFORMATEX format;

    assert(config != NULL);

    /* Multiple initialization is not possible */
    if (st_initialize_count > 0) {
        return;
    }

    /* Get config */
    st_config = (*config);

    /* Allocate buffer for decoding area */
    st_decode_buffer = (int32_t**)malloc(sizeof(int32_t*) * st_config.num_channels);
    for (i = 0; i < st_config.num_channels; i++) {
        st_decode_buffer[i] = (int32_t*)malloc(sizeof(int32_t) * DECODE_BUFFER_NUM_SAMPLES);
        memset(st_decode_buffer[i], 0, sizeof(int32_t) * DECODE_BUFFER_NUM_SAMPLES);
    }

    /* Initialize COM */
    hr = CoInitializeEx(NULL, COINIT_SPEED_OVER_MEMORY);
    assert(SUCCEEDED(hr));

    /* Get multimedia device enumerator */
    hr = CoCreateInstance(&st_CLSID_MMDeviceEnumerator, NULL, CLSCTX_ALL, &st_IID_IMMDeviceEnumerator, &device_enumerator);
    assert(SUCCEEDED(hr));

    /* Get default device */
    hr = IMMDeviceEnumerator_GetDefaultAudioEndpoint(device_enumerator, eRender, eConsole, &audio_device);
    assert(SUCCEEDED(hr));
    IMMDeviceEnumerator_Release(device_enumerator);

    /* Get client */
    hr = IMMDevice_Activate(audio_device, &st_IID_IAudioClient, CLSCTX_ALL, NULL, &audio_client);
    assert(SUCCEEDED(hr));
    IMMDevice_Release(audio_device);

    /* Output format specification */
    ZeroMemory(&format, sizeof(WAVEFORMATEX));
    format.wFormatTag = WAVE_FORMAT_PCM;
    format.nChannels = st_config.num_channels;
    format.nSamplesPerSec = st_config.sampling_rate;
    format.wBitsPerSample = st_config.bits_per_sample;
    format.nBlockAlign = (format.nChannels * format.wBitsPerSample) / 8;
    format.nAvgBytesPerSec = format.nSamplesPerSec * format.nBlockAlign;

    /* Check if the output format is supported */
    {
        WAVEFORMATEX closest_format, *pformat;
        pformat = &closest_format;
        hr = IAudioClient_IsFormatSupported(audio_client, AUDCLNT_SHAREMODE_SHARED, &format, &pformat);
        if (hr != S_OK) {
            fprintf(stderr, "Unsupported format for WASAPI Playback. \n");
            exit(1);
        }
    }

    /* Client initialization */
    hr = IAudioClient_Initialize(audio_client,
            AUDCLNT_SHAREMODE_SHARED, /* Shared mode */
            AUDCLNT_STREAMFLAGS_RATEADJUST /* Use rate conversion (play at a rate that matches the input waveform) */
            | AUDCLNT_STREAMFLAGS_AUTOCONVERTPCM /* Enable automatic insertion of rate conversions */
            | AUDCLNT_STREAMFLAGS_SRC_DEFAULT_QUALITY, /* Use good quality rate conversion */
            REQUESTED_SOUND_BUFFER_DURATION, 0, &format, NULL);
    if (FAILED(hr)) {
        fprintf(stderr, "Failed to initialize WASAPI client. \n");
        exit(2);
    }

    /* Sampling rate conversion settings */
    {
        IAudioClockAdjustment* clock_adj;
        hr = IAudioClient_GetService(audio_client, &st_IID_IAudioClockAdjustment, &clock_adj);
        assert(SUCCEEDED(hr));

        hr = IAudioClockAdjustment_SetSampleRate(clock_adj, st_config.sampling_rate);
        assert(SUCCEEDED(hr));
        IAudioClockAdjustment_Release(clock_adj);
    }

    /* Get the renderer to write the buffer to */
    hr = IAudioClient_GetService(audio_client, &st_IID_IAudioRenderClient, &audio_render_client);
    assert(SUCCEEDED(hr));

    /* Get the buffer size for writing */
    hr = IAudioClient_GetBufferSize(audio_client, &buffer_frame_size);
    assert(SUCCEEDED(hr));

    /* Start playing */
    hr = IAudioClient_Start(audio_client);
    assert(SUCCEEDED(hr));

    st_initialize_count++;

    while (1) {
        int16_t* buffer;
        /* Latency: too small and you'll get choppy results, too large and you'll get a lot of delays */
        const uint32_t buffer_latency = buffer_frame_size / 50;
        uint32_t padding_size, available_buffer_frame_size;

        /* Get padding frame size (amount of data in the sound buffer that has not yet been output) */
        hr = IAudioClient_GetCurrentPadding(audio_client, &padding_size);
        assert(SUCCEEDED(hr));

        /* Get the writable frame size */
        available_buffer_frame_size = buffer_latency - padding_size;

        /* Get buffer for writing */
        hr = IAudioRenderClient_GetBuffer(audio_render_client, available_buffer_frame_size, &buffer);
        assert(SUCCEEDED(hr));

        /* Write with interleaving. A group of samples for the number of channels is one frame */
        for (i = 0; i < available_buffer_frame_size; i++) {
            uint32_t ch;
            /* If the buffer is full, request the next data immediately */
            if (st_buffer_pos >= DECODE_BUFFER_NUM_SAMPLES) {
                st_config.sample_request_callback(st_decode_buffer, st_config.num_channels, DECODE_BUFFER_NUM_SAMPLES);
                st_buffer_pos = 0;
            }

            /* Fill the interleaved buffer with data */
            for (ch = 0; ch < st_config.num_channels; ch++) {
                *buffer++ = (int16_t)st_decode_buffer[ch][st_buffer_pos];
            }
            st_buffer_pos++;
        }

        /* Free the buffer */
        hr = IAudioRenderClient_ReleaseBuffer(audio_render_client, available_buffer_frame_size, 0);
        assert(SUCCEEDED(hr));
    }
}

/* End. Release the initialized resources here. */
void SRLAPlayer_Finalize(void)
{
    if (st_initialize_count == 1) {
        IAudioClient_Stop(audio_client);
        IAudioClient_Release(audio_client);
        IAudioRenderClient_Release(audio_render_client);
    }

    st_initialize_count--;
}
