#ifndef SRLAPLAYER_H_INCLUDED
#define SRLAPLAYER_H_INCLUDED

#include <stdint.h>

/* Output request callback */
typedef void (*SRLASampleRequestCallback)(
        int32_t **buffer, uint32_t num_channels, uint32_t num_samples);

/* Player initialization config */
struct SRLAPlayerConfig {
    uint32_t sampling_rate;
    uint16_t num_channels;
    uint16_t bits_per_sample;
    SRLASampleRequestCallback sample_request_callback;
};

#ifdef __cplusplus
extern "C" {
#endif

/* Initialization This function initializes the device driver and starts playback. */
void SRLAPlayer_Initialize(const struct SRLAPlayerConfig *config);

/* End. Release the initialized resources here. */
void SRLAPlayer_Finalize(void);

#ifdef __cplusplus
}
#endif

#endif /* SRLAPLAYER_H_INCLUDED */
