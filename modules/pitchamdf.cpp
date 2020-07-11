#include <math.h>
#include <cstdint>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include "pitchamdf.h"

#define SPFLOAT float
#define SP_OK 1
#define SP_NOT_OK 0

using namespace daisysp;

int sp_auxdata_alloc(void **aux, size_t size)
{
    *aux = malloc(size);
    memset(*aux, 0, size);
    return SP_OK;
}

int PitchAmdf::Init(float sample_rate, float min_freq, float max_freq)
{
    SPFLOAT srate, downs;
    int32_t size, minperi, maxperi, downsamp, upsamp, msize, bufsize;
    uint32_t interval;

    this->imincps = min_freq;
    this->imaxcps = max_freq;

    /* TODO: should we expose these variables? */
    this->icps = 0;
    this->imedi = 1;
    this->idowns = 1;
    this->iexcps = 0;
    this->irmsmedi = 0;

    this->inerr = 0;
    downs = this->idowns;

    if (downs < (-1.9)) {
        upsamp = (int)lrintf((-downs));
        downsamp = 0;
        srate = sample_rate * (SPFLOAT)upsamp;
    } else {
        downsamp = (int)lrintf(downs);
        if (downsamp < 1) downsamp = 1;
        srate = sample_rate / (SPFLOAT)downsamp;
        upsamp = 0;
    }

    minperi = (int32_t)(srate / this->imaxcps);
    maxperi = (int32_t)(0.5 + srate / this->imincps);
    if (maxperi <= minperi) {
        this->inerr = 1;
        return SP_NOT_OK;
    }

    if (this->iexcps < 1)
        interval = maxperi;
    else
        interval = (uint32_t)(srate / this->iexcps);

    size = maxperi + interval;
    bufsize = sizeof(SPFLOAT)*(size + maxperi + 2);

    this->srate = srate;
    this->downsamp = downsamp;
    this->upsamp = upsamp;
    this->minperi = minperi;
    this->maxperi = maxperi;
    this->size = size;
    this->readp = 0;
    this->index = 0;
    this->lastval = 0.0;

    if (this->icps < 1) {
        this->peri = (minperi + maxperi) / 2;
    } else {
        this->peri = (int)(srate / this->icps);
    }

    if (this->irmsmedi < 1) {
        this->rmsmedisize = 0;
    } else {
        this->rmsmedisize = ((int)lrintf(this->irmsmedi))*2+1;
    }

    this->rmsmediptr = 0;

    if (this->rmsmedisize) {
        msize = this->rmsmedisize * 3 * sizeof(SPFLOAT);
        sp_auxdata_alloc(&this->rmsmedian, msize);
    }

    if (this->imedi < 1) {
        this->medisize = 0;
    } else {
        this->medisize = (int)lrintf(this->imedi) * 2 + 1;
    }

    this->mediptr = 0;

    if (this->medisize) {
        msize = this->medisize * 3 * sizeof(SPFLOAT);
        sp_auxdata_alloc(&this->median, msize);
    }

    sp_auxdata_alloc(&this->buffer, bufsize);
    return SP_OK;
}

#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp

static SPFLOAT medianvalue(uint32_t n, SPFLOAT *vals)
{   
    /* vals must point to 1 below relevant data! */
    uint32_t i, ir, j, l, mid;
    uint32_t k = (n + 1) / 2;
    SPFLOAT a, temp;

    l = 1;
    ir = n;
    while (1) {
        if (ir <= l+1) {
            if (ir == l+1 && vals[ir] < vals[l]) {
                SWAP(vals[l], vals[ir]);
            }
            return vals[k];
        } else {
            mid = (l+ir) >> 1;
            SWAP(vals[mid], vals[l+1]);
            if (vals[l+1] > vals[ir]) {
                SWAP(vals[l+1], vals[ir]);
            }
            if (vals[l] > vals[ir]) {
                SWAP(vals[l], vals[ir]);
            }
            if (vals[l+1] > vals[l]) {
                SWAP(vals[l+1], vals[l]);
            }
            i = l + 1;
            j = ir;
            a = vals[l];
            while (1) {
                do i++; while (vals[i] < a);
                do j--; while (vals[j] > a);
                if (j < i) break;
                SWAP(vals[i], vals[j]);
            }
            vals[l] = vals[j];
            vals[j] = a;
            if (j >= k) ir = j-1;
            if (j <= k) l = i;
        }
    }
}
#undef SWAP

float PitchAmdf::Process(float sample)
{
    SPFLOAT *buffer = (SPFLOAT*)this->buffer;
    SPFLOAT *rmsmedian = (SPFLOAT*)this->rmsmedian;
    int32_t rmsmedisize = this->rmsmedisize;
    int32_t rmsmediptr = this->rmsmediptr;
    SPFLOAT *median = (SPFLOAT*)this->median;
    int32_t medisize = this->medisize;
    int32_t mediptr = this->mediptr;
    int32_t size = this->size;
    int32_t index = this->index;
    int32_t minperi = this->minperi;
    int32_t maxperi = this->maxperi;
    SPFLOAT srate = this->srate;
    int32_t peri = this->peri;
    int32_t upsamp = this->upsamp;
    SPFLOAT upsmp = (SPFLOAT)upsamp;
    SPFLOAT lastval = this->lastval;
    SPFLOAT newval, delta;
    int32_t readp = this->readp;
    int32_t interval = size - maxperi;
    int i;
    int32_t i1, i2;
    SPFLOAT val, rms;
    SPFLOAT sum;
    SPFLOAT acc, accmin, diff;

    if (upsamp) {
        newval = sample;
        delta = (newval-lastval) / upsmp;
        lastval = newval;

        for (i=0; i<upsamp; i++) {
            newval += delta;
            buffer[index++] = newval;

            if (index == size) {
                peri = minperi;
                accmin = 0.0;
                for (i2 = 0; i2 < size; ++i2) {
                    diff = buffer[i2+minperi] - buffer[i2];
                    if (diff > 0) accmin += diff;
                    else accmin -= diff;
                }
                for (i1 = minperi + 1; i1 <= maxperi; ++i1) {
                    acc = 0.0;
                    for (i2 = 0; i2 < size; ++i2) {
                        diff = buffer[i1+i2] - buffer[i2];
                        if (diff > 0) acc += diff;
                        else acc -= diff;
                        if (acc > accmin) {
                            break;
                        }
                    }
                    if (acc < accmin) {
                        accmin = acc;
                        peri = i1;
                    }
                }

                for (i1 = 0; i1 < interval; i1++) { 
                    buffer[i1] = buffer[i1+interval]; 
                }

                index = maxperi;

                if (medisize) {
                    median[mediptr] = (SPFLOAT)peri;
                    for (i1 = 0; i1 < medisize; i1++) {
                        median[medisize+i1] = median[i1];
                    }

                    median[medisize*2+mediptr] =
                    medianvalue(medisize, &median[medisize-1]);
                    peri = (int32_t)median[medisize*2 +
                        ((mediptr+medisize/2+1) % medisize)];

                    mediptr = (mediptr + 1) % medisize;
                    this->mediptr = mediptr;
                }
            }
        }
        this->lastval = lastval;
    } else {
        int32_t  downsamp = this->downsamp;
        buffer[index++] = sample;
        readp += downsamp;

        if (index == size) {
            peri = minperi;
            accmin = 0.0;

            for (i2 = 0; i2 < size; ++i2) {
                diff = buffer[i2+minperi] - buffer[i2];
                if (diff > 0.0) accmin += diff;
                else accmin -= diff;
            }

            for (i1 = minperi + 1; i1 <= maxperi; ++i1) {
                acc = 0.0;
                for (i2 = 0; i2 < size; ++i2) {
                    diff = buffer[i1+i2] - buffer[i2];
                    if (diff > 0.0) acc += diff;
                    else acc -= diff;
                    if (acc > accmin) {
                        break;
                    }
                }
                if (acc < accmin) {
                    accmin = acc;
                    peri = i1;
                }
            }

            for (i1 = 0; i1 < interval; i1++) {
                buffer[i1] = buffer[i1+interval];
            }

            index = maxperi;

            if (medisize) {
                median[mediptr] = (SPFLOAT)peri;

                for (i1 = 0; i1 < medisize; i1++) {
                    median[medisize+i1] = median[i1];
                }

                median[medisize*2+mediptr] =
                medianvalue(medisize, &median[medisize-1]);
                peri = (int32_t)median[medisize*2 +
                    ((mediptr+medisize/2+1) % medisize)];

                mediptr = (mediptr + 1) % medisize;
                this->mediptr = mediptr;
            }
        }
    }
    buffer = &buffer[(index + size - peri) % size];
    sum = 0.0;
    for (i1=0; i1<peri; i1++) {
        val = buffer[i1];
        sum += (SPFLOAT)(val * val);
    }
    if (peri==0)      
        rms = 0.0;
    else
        rms = (SPFLOAT)sqrt(sum / (SPFLOAT)peri);
    if (rmsmedisize) {
        rmsmedian[rmsmediptr] = rms;
        for (i1 = 0; i1 < rmsmedisize; i1++) {
            rmsmedian[rmsmedisize+i1] = rmsmedian[i1];
        }

        rmsmedian[rmsmedisize*2+rmsmediptr] =
            medianvalue(rmsmedisize, &rmsmedian[rmsmedisize-1]);
        rms = rmsmedian[rmsmedisize*2 +
            ((rmsmediptr+rmsmedisize/2+1) % rmsmedisize)];

        rmsmediptr = (rmsmediptr + 1) % rmsmedisize;
        this->rmsmediptr = rmsmediptr;
    }

    if (peri==0) {
        cps_ = 0.0;
    } else {
        cps_ = srate / (SPFLOAT)peri;
    }

    rms_ = rms;
    this->index = index;
    this->peri = peri;
    this->readp = readp;

    return cps_;
}

// int main() {
//     PitchAmdf pitchAmdf;
//     pitchAmdf.Init(48000, 200.f, 500.f);

//     SPFLOAT frequency = 350.f;
//     SPFLOAT value = 0.f;
//     SPFLOAT output = 0.f;
//     for(int i=0;i<48000;i++) {
//         value = sin((double)i / 48000 * M_PI * 2.0 * frequency);
//         output = pitchAmdf.Process(value);
//     }

//     printf("done, %f %f\n", output, pitchAmdf.GetRms());

//     return 0;
// }