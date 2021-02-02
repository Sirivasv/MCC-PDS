[ydoa1, fs] = audioread("snap_out_file_doa_1.wav")
[ydoa2, fs] = audioread("snap_out_file_doa_2.wav")
[ypris1, fs] = audioread("snap_pristine_clean_channel_1.wav")
[ypris2, fs] = audioread("snap_pristine_clean_channel_2.wav")
%[ydoa1, fs] = audioread("snap_pristine_clean_channel_1.wav")
%[ydoa2, fs] = audioread("snap_pristine_clean_channel_2.wav")

[nc1, nr1] = size(ydoa1)
[nc2, nr1] = size(ydoa2)
[nc3, nr1] = size(ypris1)
[nc4, nr1] = size(ypris2)
[mincn, imincn] = min([nc1,nc2,nc3,nc4])

N = mincn;

ydoa1 = ydoa1(1:N)
ydoa1 = reshape(ydoa1, 1,N)

ydoa2 = ydoa2(1:N)
ydoa2 = reshape(ydoa2, 1,N)

ypris1 = ypris1(1:N)
ypris1 = reshape(ypris1, 1,N)

ypris2 = ypris2(1:N)
ypris2 = reshape(ypris2, 1,N)

[SDR,SIR,SAR,perm]=bss_eval_sources([ydoa1;ydoa2],[ypris1;ypris2])