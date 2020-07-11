#pragma once
#ifndef DSY_PITCH_AMDF_H
#define DSY_PITCH_AMDF_H
#ifdef __cplusplus

#define SPFLOAT float

namespace daisysp
{
/** placeholder: Generates a normalized signal moving from 0-1 at the specified frequency.

\todo placeholder: Selecting which channels should be initialized/included in the sequence conversion.
\todo placeholder: Setup a similar start function for an external mux, but that seems outside the scope of this file.

*/
class PitchAmdf
{
  public:
    PitchAmdf() {}
    ~PitchAmdf() {}
    /** Initializes the PitchAmdf module
    sample rate, minFreq, maxFreq are in Hz
    Additional Init functions have defaults when arg is not specified:
    - minFreq = 200.0f
    - maxFreq = 500.0f
    */
    int Init(float sample_rate, float min_freq, float max_freq);

    /** Initialize PitchAmdf with samplerate and freq
    */
    inline int Init(float sample_rate, float min_freq)
    {
      return Init(sample_rate, min_freq, 500.0f);
    }

    /** Initialize PitchAmdf with samplerate
    */
    inline int Init(float sample_rate) { 
      return Init(sample_rate, 200.0f);
    }
    /** processes PitchAmdf and returns detected frequency
    */
    float Process(float sample);

    /** get volume of detected pitch
    */
    float GetRms() {
      return rms_;
    }

    float GetCps() {
      return cps_;
    }

  private:
    float rms_, cps_;
    // ported from soundpipe
    SPFLOAT imincps, imaxcps, icps;
    SPFLOAT imedi, idowns, iexcps, irmsmedi;
    SPFLOAT srate;
    SPFLOAT lastval;
    int32_t downsamp;
    int32_t upsamp;
    int32_t minperi;
    int32_t maxperi;
    int32_t index;
    int32_t readp;
    int32_t size;
    int32_t peri;
    int32_t medisize;
    int32_t mediptr;
    int32_t rmsmedisize;
    int32_t rmsmediptr;
    int inerr;
    void *median;
    void *rmsmedian;
    void *buffer;
};
} // namespace daisysp
#endif
#endif
