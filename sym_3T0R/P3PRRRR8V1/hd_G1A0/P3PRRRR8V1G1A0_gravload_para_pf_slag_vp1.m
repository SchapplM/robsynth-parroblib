% Calculate Gravitation load for parallel robot
% P3PRRRR8V1G1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4,alpha2,d2,d4,theta1]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [3x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 16:50
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3PRRRR8V1G1A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V1G1A0_gravload_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V1G1A0_gravload_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3PRRRR8V1G1A0_gravload_para_pf_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V1G1A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR8V1G1A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRRR8V1G1A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V1G1A0_gravload_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V1G1A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 16:49:48
% EndTime: 2020-08-06 16:49:49
% DurationCPUTime: 1.07s
% Computational Cost: add. (531->132), mult. (1289->259), div. (36->7), fcn. (1074->22), ass. (0->116)
t926 = cos(pkin(3));
t927 = legFrame(3,3);
t913 = sin(t927);
t916 = cos(t927);
t897 = -t913 * g(1) + t916 * g(2);
t900 = t916 * g(1) + t913 * g(2);
t923 = sin(pkin(6));
t925 = cos(pkin(6));
t959 = t897 * t925 - t900 * t923;
t924 = sin(pkin(3));
t993 = g(3) * t924;
t1002 = t959 * t926 + t993;
t931 = sin(qJ(2,3));
t937 = cos(qJ(2,3));
t958 = t897 * t923 + t900 * t925;
t947 = t1002 * t931 + t958 * t937;
t928 = legFrame(2,3);
t914 = sin(t928);
t917 = cos(t928);
t898 = -t914 * g(1) + t917 * g(2);
t901 = t917 * g(1) + t914 * g(2);
t957 = t898 * t925 - t901 * t923;
t1003 = t957 * t926 + t993;
t933 = sin(qJ(2,2));
t939 = cos(qJ(2,2));
t956 = t898 * t923 + t901 * t925;
t946 = t1003 * t933 + t956 * t939;
t929 = legFrame(1,3);
t915 = sin(t929);
t918 = cos(t929);
t899 = -t915 * g(1) + t918 * g(2);
t902 = t918 * g(1) + t915 * g(2);
t955 = t899 * t925 - t902 * t923;
t1004 = t955 * t926 + t993;
t935 = sin(qJ(2,1));
t941 = cos(qJ(2,1));
t954 = t899 * t923 + t902 * t925;
t945 = t1004 * t935 + t954 * t941;
t992 = g(3) * t926;
t936 = cos(qJ(3,3));
t930 = sin(qJ(3,3));
t1000 = pkin(2) * t930;
t997 = pkin(2) * t936;
t903 = -t937 * pkin(5) + t931 * t997;
t965 = 0.1e1 / (t926 * t1000 + t903 * t924);
t1007 = 0.1e1 / t936 * t965;
t938 = cos(qJ(3,2));
t996 = pkin(2) * t938;
t904 = -t939 * pkin(5) + t933 * t996;
t932 = sin(qJ(3,2));
t999 = pkin(2) * t932;
t964 = 0.1e1 / (t904 * t924 + t926 * t999);
t1006 = 0.1e1 / t938 * t964;
t940 = cos(qJ(3,1));
t995 = pkin(2) * t940;
t905 = -t941 * pkin(5) + t935 * t995;
t934 = sin(qJ(3,1));
t998 = pkin(2) * t934;
t963 = 0.1e1 / (t905 * t924 + t926 * t998);
t1005 = 0.1e1 / t940 * t963;
t1001 = m(3) / pkin(2);
t919 = m(1) + m(2) + m(3);
t994 = g(3) * t919;
t991 = rSges(3,2) * t924;
t990 = t924 * t936;
t989 = t924 * t938;
t988 = t924 * t940;
t987 = t926 * t931;
t986 = t926 * t933;
t985 = t926 * t935;
t882 = -t923 * t913 + t916 * t925;
t984 = t931 * t882;
t885 = t925 * t913 + t916 * t923;
t983 = t931 * t885;
t883 = -t923 * t914 + t917 * t925;
t982 = t933 * t883;
t886 = t925 * t914 + t917 * t923;
t981 = t933 * t886;
t884 = -t923 * t915 + t918 * t925;
t980 = t935 * t884;
t887 = t925 * t915 + t918 * t923;
t979 = t935 * t887;
t978 = t937 * t882;
t977 = t937 * t885;
t976 = t939 * t883;
t975 = t939 * t886;
t974 = t941 * t884;
t973 = t941 * t887;
t972 = rSges(3,2) * t992;
t971 = (((t959 * t924 - t992) * rSges(3,1) + t947 * rSges(3,2)) * t936 + (t947 * rSges(3,1) - t959 * t991 + t972) * t930) * t1007;
t970 = (((t957 * t924 - t992) * rSges(3,1) + t946 * rSges(3,2)) * t938 + (t946 * rSges(3,1) - t957 * t991 + t972) * t932) * t1006;
t969 = (((t955 * t924 - t992) * rSges(3,1) + t945 * rSges(3,2)) * t940 + (t945 * rSges(3,1) - t955 * t991 + t972) * t934) * t1005;
t912 = m(2) * rSges(2,2) - m(3) * rSges(3,3);
t943 = m(2) * rSges(2,1);
t968 = (t947 * t912 + (-t1002 * t937 + t931 * t958) * (t943 + (rSges(3,1) * t936 - rSges(3,2) * t930) * m(3))) * t1007;
t967 = (t946 * t912 + (-t1003 * t939 + t933 * t956) * (t943 + (rSges(3,1) * t938 - rSges(3,2) * t932) * m(3))) * t1006;
t966 = (t945 * t912 + (-t1004 * t941 + t935 * t954) * (t943 + (rSges(3,1) * t940 - rSges(3,2) * t934) * m(3))) * t1005;
t953 = t924 * t1000 - t903 * t926;
t952 = -t904 * t926 + t924 * t999;
t951 = -t905 * t926 + t924 * t998;
t908 = pkin(5) * t935 + t941 * t995;
t907 = pkin(5) * t933 + t939 * t996;
t906 = pkin(5) * t931 + t937 * t997;
t893 = -t923 * t985 + t925 * t941;
t892 = -t923 * t986 + t925 * t939;
t891 = -t923 * t987 + t925 * t937;
t890 = t923 * t941 + t925 * t985;
t889 = t923 * t939 + t925 * t986;
t888 = t923 * t937 + t925 * t987;
t872 = t923 * t908 - t951 * t925;
t871 = t923 * t907 - t952 * t925;
t870 = t923 * t906 - t953 * t925;
t869 = t908 * t925 + t951 * t923;
t868 = t907 * t925 + t952 * t923;
t867 = t906 * t925 + t953 * t923;
t1 = [((-t890 * t918 - t915 * t893) * t934 - t884 * t988) * t966 + ((-t889 * t917 - t914 * t892) * t932 - t883 * t989) * t967 + ((-t888 * t916 - t913 * t891) * t930 - t882 * t990) * t968 - m(4) * g(1) + (-(t869 * t918 - t915 * t872) * t963 - (t868 * t917 - t914 * t871) * t964 - (t867 * t916 - t913 * t870) * t965) * t994 + ((-(t926 * t974 - t979) * t995 - pkin(5) * (t926 * t980 + t973)) * t969 + (-(t926 * t976 - t981) * t996 - pkin(5) * (t926 * t982 + t975)) * t970 + (-(t926 * t978 - t983) * t997 - pkin(5) * (t926 * t984 + t977)) * t971) * t1001; ((-t915 * t890 + t893 * t918) * t934 - t887 * t988) * t966 + ((-t914 * t889 + t892 * t917) * t932 - t886 * t989) * t967 + ((-t913 * t888 + t891 * t916) * t930 - t885 * t990) * t968 - m(4) * g(2) + (-(t869 * t915 + t872 * t918) * t963 - (t868 * t914 + t871 * t917) * t964 - (t867 * t913 + t870 * t916) * t965) * t994 + ((-(t926 * t973 + t980) * t995 - (t926 * t979 - t974) * pkin(5)) * t969 + (-(t926 * t975 + t982) * t996 - (t926 * t981 - t976) * pkin(5)) * t970 + (-(t926 * t977 + t984) * t997 - (t926 * t983 - t978) * pkin(5)) * t971) * t1001; (-m(4) - 0.3e1 * t919) * g(3);];
taugX  = t1;
