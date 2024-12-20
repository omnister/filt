# filt
linux command line c program to digitally filter a stream of data in time domain

The program synthesizes digital filters as cascaded biquads and
interprets them to produce low-pass, high-pass, band-pass and
band-reject filters to operate on ascii data streams. 

Supported filter types are: Butterworth, maxflat, chebyshev, modified
chebyshev, inverse chebyshev, and bessel. 


-------------------------- MAN PAGE ------------------------------------------
lpf(l)                                                                  lpf(l)

NAME
       lpf, hpf, bpf, brf - filter a stream of input data in the time domain

SYNOPSIS
       lpf [ -dLqxX?  ] [ -f filename ] [ -t filter_type ] [ -n filter_order ]
       [ -c center_frequency ] [ -b filter_bandwidth ] [ -r ripple_tolerance ]
       [ -s shape_factor ] [ -F fnorm_char ] [ -i dc_value ] [ -I nsamples ] [
       -S nsamples ]

DESCRIPTION
       All of the Filter suite of programs read "time,value" pairs from  stan‐
       dard  input and write a filtered set of values to standard output.  lpf
       is a lowpass filter, hpf is a highpass filter, bpf is a  bandpass  fil‐
       ter,  and  brf is a bandreject filter.  All of the programs assume that
       the input data is evenly spaced in time. In fact, the "time"  value  is
       ignored  and  merely passed unchanged to the output.  The output format
       is ASCII "time,value" pairs, with one pair per line,  and  is  suitable
       for  input  to pdplot program.  Lines which do not begin with a numeric
       value are passed directly to the output, so  that  formatting  commands
       can be included in the data.  The following options are recognized:

       -t  filter_type
               Filter_type  can  be  either  "butterworth", "maxflat", "cheby‐
              shev", "modcheb", "invcheb", or "bessel".  Default  is  "butter‐
              worth".   The Butterworth approximation is a special case of the
              maximally-flat approximation with its loss at the  edge  of  the
              passband  set to 3 dB.  The modified Chebyshev approximation ap‐
              plies to even-order filters only, shifting their lowest  reflec‐
              tion zeros (in the LP case) to zero.  Hence, even-order modified
              Chebyshev filters have unity gain at dc.  Inverse Chebyshev fil‐
              ters have a maximally-flat passband and equiminima stopband with
              finite transmission zeros.  They have less total delay and delay
              distortion  than  a  comparable Chebyshev filter.  They are also
              particularly useful in this program for bandreject filters.

       -n  filter_order
              This gives the number of natural  frequencies  in  the  low-pass
              prototype.   Bandpass  and  bandreject filters will actually use
              twice the specified number of frequencies.  Default  is  8th-or‐
              der.

       -c  center_frequency
              Specifies  the  the  arithmetic  center frequency of bandpass or
              bandreject filters  relative  to  the  Nyquist  frequency  (sam‐
              pling_frequency  /  2).  The parameter is ignored for lpf or hpf
              filters.  Default is 0.5.

       -b  bandwidth
              Specifies the filter bandwidth in normalized frequency  relative
              to  Nyquist.  For example, if the sampling rate is 1 kHz and the
              desired low-pass cutoff frequency is 50 Hz, the bandwidth should
              be specified as "-b 0.1".  Default is 0.1.  Note: Bessel filters
              are not specified in terms of bandwidth, but in  terms  of  gain
              and  phase  flatness  for  a given delay.  The number of samples
              that a Bessel filter delays the data  is  equal  to  1/(PI*band‐
              width),  as  specified on the command line.  Higher order Bessel
              filters do not increase the delay, they just produce a more lin‐
              ear delay characteristic out to higher frequencies.

       -r  ripple_tolerance
              Specify  the amount of passband ripple, in dB, for the Chebyshev
              approximation.  Also specifies the loss at the edge of the pass‐
              band for maximally-flat and inverse Chebyshev filters.  This pa‐
              rameter is ignored for Butterworth filters.  Default  is  3.0103
              dB.

       -s  shape_factor
              This  parameter  applies  to inverse Chebyshev filters only, and
              specifies the ratio of the stopband edge frequency to the  pass‐
              band edge frequency (Fsb/Fpb) for lowpass filters.  It specifies
              Fpb/Fsb for highpass filters, delta_Fsb/delta_Fpb  for  bandpass
              filters,  and  delta_Fpb/delta_Fsb  for bandreject filters.  De‐
              fault  shape_factor is 2.

       -F  fnorm_char
              Select the frequency normalization factor.  This affects the in‐
              terpretation  of  the -cbxXLP options.  If  fnorm_char is 'n' or
              'N', then normalize by the Nyquist frequency (sampling_frequency
              /  2).   If  fnorm_char is 's' or 'S', normalize by the sampling
              frequency.  Default is to normalize by the Nyquist frequency.

       -i  dc_value
              Specify the filter's initial input value, i.e., the filter  acts
              as if the filter input had been  dc_value from t = (minus infin‐
              ity) to t = 0.  Of course, the default value is 0.

       -x     Causes the program to output the calculated  frequency  response
              in  pdplot(1)  format.   No  filtering  is done.  This option is
              helpful in choosing the filter parameters before actually  doing
              the time domain filtering.

       -X     Like the -x option, with the filter's delay plotted as well.

       -L     Like the -X option, plotting the more traditional loss and delay
              curves.

       -P     A variant on the -L option, wherein only the passband  is  plot‐
              ted.

       -I  nsamples
              Outputs,  in  pdplot(1) format,  nsamples of the time-domain im‐
              pulse response of the filter.  This and the -S option are useful
              to  observe the settling behavior of the filter and to make sure
              that the program has correctly constructed a  stable filter :-).

       -S  nsamples
              Outputs, in pdplot(1) format,  nsamples of the time-domain  step
              response of the filter.

       -d     Causes the program to dump the calculated z-domain biquad param‐
              eters to standard output.  These parameters can be used  to  du‐
              plicate  the filter in hardware.  (Note, however, that this pro‐
              gram's scaling, pairing, and ordering of pole-zero pairs is  not
              optimized  for  fixed-point arithmetic.)  They can also be saved
              and read back in with the -f option.  No filtering is done  when
              this option is specified.

       -f  filename
              Causes  the  program to load the z-domain biquad parameters from
              filename and to use the  parameters  for  subsequent  filtering.
              This  option  allows parameters to be imported from other filter
              design programs.  Also, two or more dump files created with  the
              -d  option  can be concatenated to form composite filters.  When
              using imported parameters, the -tncbrs options are all ignored.

       -q     do not pass the non-data lines to the output.

       -?     print out a quick summary of options and exit.

EXAMPLES
       If the sampling rate is 1 kHz and the desired 3 dB frequency is 50  Hz,
       a 20th-order Butterworth lowpass filter would be invoked as:

           "lpf -b0.1 -n20 -tbutt < input_file > output_file"

       To  see the frequency domain performance of this filter, use the -X op‐
       tion:

           "lpf -b0.1 -n20 -tbutt -X | pdplot"

       This example compares the magnitude response of  an  inverse  Chebyshev
       filter to a commensurate plain Chebyshev filter:

           "(echo xset 0.4 0.6; bpf -tcheb -x; bpf -tinvcheb -s2 -x) | pdplot"

       The phase distortion of the Bessel approximation (see BUGS) is revealed
       by this example:

           "(lpf -tbess -b.1 -n10  -S100;  lpf  -tbess  -b.01  -n10  -S100)  |
       pdplot"

       Linear phase filtering can be done by running the data through the same
       filter in both forward and reverse time.  This can be done easily  with
       the tac(1) program that "cats" lines in reverse order:

           "tac datafile | lpf > temp"
           "tac temp | lpf > outputfile"

       Or try this, to verify the symmetry of the filter's impulse response:

           nawk 'BEGIN {
             for (i=199; i>=0; i--)
               if (i == 100) print i, 1; else print i, 0
           }' | lpf | tac | lpf | pdplot

       This combines forward/backward filtering with the -i option to obtain a
       symmetric step:

         nawk 'BEGIN {
           for (i=199; i>=0; i--)
             if (i >= 100) print i, 1; else print i, 0
         }' | lpf -i1 | tac | lpf | pdplot

       The -i option is illustrated in the following example, where  the  step
       response  transient  is suppressed by a priori knowledge of the initial
       condition:

           "(lpf -S100; lpf -i1.0 -S100) | pdplot"

       An interesting DSP principle is demonstrated in this example:

           "bpf -tinvcheb -n5 -b.1 -I80 | pdplot"

       A symmetric bandpass filter, centered on half  the  Nyquist  frequency,
       has  an impulse response where every other sample is zero.  This can be
       exploited to simplify its hardware implementation.

SEE ALSO
       pdplot(1)

AUTHOR
       Written by Rick Walker, HPL.

       Modifications by Scott Willingham, HPL.

BUGS
       The Bessel filter option does not create a filter  with  linear  phase!
       This  results  from  the bilinear transformation of the s-domain Bessel
       polynomial into the the z-domain.  The frequency mapping is  inherently
       warped,  distorting  the  phase and delay response.  Narrowband lowpass
       Bessel filters (i.e., with long delay) should be acceptable.  Highpass,
       bandpass,  and  bandreject  transformations all destroy linear phase as
       well.  The program allows such filters to be specified,  but  issues  a
       warning on stderr.
