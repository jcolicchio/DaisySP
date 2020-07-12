#pragma once
#ifndef DSY_WAVEFORMBUFFER_H
#define DSY_WAVEFORMBUFFER_H

#include <stdint.h>
#ifdef __cplusplus

#define WAVEFORM_MAX_BUFFER_LENGTH (480)
#define WAVEFORM_SAMPLE_MIN (-1.0)
#define WAVEFORM_SAMPLE_MAX (1.0)
#define WAVEFORM_BUFFER_START (2)

namespace daisysp
{
/** Buffer a waveform for visualization
    @author Joseph Colicchio
    @date   2020
*/

class WaveformBuffer {
public:
    WaveformBuffer() {
        _isActive = false;
    }
    ~WaveformBuffer() {}
    void Init(uint32_t sampleRate, double frequency, uint32_t refreshRate, uint32_t maxBufferLength) {
        Init(sampleRate, frequency, refreshRate, maxBufferLength, WAVEFORM_SAMPLE_MIN, WAVEFORM_SAMPLE_MAX);
    }
    void Init(uint32_t sampleRate, double frequency, uint32_t refreshRate, uint32_t maxBufferLength, float sampleMin, float sampleMax) {
        _isActive = true;
        _sampleRate = sampleRate;
        _frequency = frequency;
        _refreshRate = refreshRate;
        _maxBufferLength = std::min(maxBufferLength, (uint32_t)WAVEFORM_MAX_BUFFER_LENGTH);
        setSampleRange(sampleMin, sampleMax);
        Trigger();
    }
    // manually start buffering, resetting everything and recording samples immediately
    void Trigger() {
        _isWriting = true;
        _time = 0.0;
        _refreshTime = 0.0;
        _fidelity = 0.0;
        _bufferIndex = WAVEFORM_BUFFER_START;
    }
    // process every audio sample at the initialized rate
    // if Process returns true, the buffer is full and it's time to send buffer() with bufferSize() length
    bool Process(float sample) {
        if (!_isActive) {
            return false;
        }
        
        if (!_isReadyToWrite() && !_isWriting) {
            return false;
        }

        _isWriting = true;
        _fidelity += 1.0;
        if (_fidelity < optimalFidelity()) {
            return false;
        }

        _fidelity -= optimalFidelity();
        _writeBuffer[_bufferIndex] = (sample - _sampleMin) / (_sampleMax - _sampleMin) * 255.0;
        _bufferIndex += 1;

        if (_bufferIndex < effectiveBufferLength() + WAVEFORM_BUFFER_START) {
            return false;
        }
        
        _bufferIndex = WAVEFORM_BUFFER_START;
        _isWriting = false;
        uint8_t *temp = _writeBuffer;
        _writeBuffer = _readBuffer;
        _readBuffer = temp;
        ((uint16_t *)_readBuffer)[0] = effectiveBufferLength();
        return true;
    }
    // the buffer to send if Process returns true
    uint8_t *buffer() {
        return _readBuffer;
    }
    // the length of buffer() to send if Process returns true
    size_t bufferSize() {
        return effectiveBufferLength() * sizeof(_readBuffer[0]) + WAVEFORM_BUFFER_START;
    }

    double frequency() {
        return _frequency;
    }
    void setFrequency(double frequency) {
        _frequency = frequency;
    }

    float sampleMin() {
        return _sampleMin;
    }
    void setSampleMin(float sampleMin) {
        _sampleMin = sampleMin;
    }

    float sampleMax() {
        return _sampleMax;
    }
    void setSampleMax(float sampleMax) {
        _sampleMax = sampleMax;
    }

    void setSampleRange(float sampleMin, float sampleMax) {
        _sampleMin = sampleMin;
        _sampleMax = sampleMax;
        if (_sampleMin > _sampleMax) {
            _sampleMin = WAVEFORM_SAMPLE_MIN;
            _sampleMax = WAVEFORM_SAMPLE_MAX;
        }
    }

private:
    uint8_t _buffer1[WAVEFORM_MAX_BUFFER_LENGTH];
    uint8_t _buffer2[WAVEFORM_MAX_BUFFER_LENGTH];
    uint8_t *_writeBuffer = _buffer1;
    uint8_t *_readBuffer = _buffer2;

    uint32_t _sampleRate;
    double _frequency;
    uint32_t _refreshRate;
    uint32_t _maxBufferLength;

    bool _isActive;
    bool _isWriting;
    double _time;
    double _refreshTime;
    double _fidelity;
    uint32_t _bufferIndex;
    float _sampleMin;
    float _sampleMax;

    bool _isReadyToWrite() {
        _time += (double)_frequency / (double)_sampleRate;
        _refreshTime += (double)_refreshRate / (double)_sampleRate;
        if (_time < 1.0) {
            return false;
        }
        _time -= 1.0;
        if (_refreshTime < 1.0) {
            return false;
        }
        _refreshTime -= 1.0;
        return true;
    }
    double maxSamples() {
        return (double)_sampleRate / (double)_frequency;
    }
    double optimalFidelity() {
        return std::max(maxSamples() / (double)_maxBufferLength, 1.0);
    }
    uint32_t effectiveBufferLength() {
        return maxSamples() / optimalFidelity();
    }
};
} // namespace daisysp
#endif
#endif