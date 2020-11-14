% Calculate Gravitation load for parallel robot
% P3PRRRR8V1G2A0
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
% Datum: 2020-08-06 17:03
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3PRRRR8V1G2A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V1G2A0_gravload_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V1G2A0_gravload_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3PRRRR8V1G2A0_gravload_para_pf_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V1G2A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR8V1G2A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRRR8V1G2A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V1G2A0_gravload_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V1G2A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 17:02:40
% EndTime: 2020-08-06 17:02:42
% DurationCPUTime: 1.18s
% Computational Cost: add. (675->134), mult. (1692->260), div. (54->7), fcn. (1311->22), ass. (0->116)
t935 = legFrame(3,2);
t921 = sin(t935);
t924 = cos(t935);
t902 = g(1) * t921 + g(2) * t924;
t932 = sin(pkin(3));
t934 = cos(pkin(3));
t905 = g(1) * t924 - g(2) * t921;
t931 = sin(pkin(6));
t920 = g(3) * t931;
t933 = cos(pkin(6));
t968 = t905 * t933 - t920;
t1025 = t902 * t932 + t968 * t934;
t939 = sin(qJ(2,3));
t1030 = t1025 * t939;
t936 = legFrame(2,2);
t922 = sin(t936);
t925 = cos(t936);
t903 = g(1) * t922 + g(2) * t925;
t906 = g(1) * t925 - g(2) * t922;
t966 = t906 * t933 - t920;
t1026 = t903 * t932 + t966 * t934;
t941 = sin(qJ(2,2));
t1029 = t1026 * t941;
t937 = legFrame(1,2);
t923 = sin(t937);
t926 = cos(t937);
t904 = g(1) * t923 + g(2) * t926;
t907 = g(1) * t926 - g(2) * t923;
t964 = t907 * t933 - t920;
t1027 = t904 * t932 + t964 * t934;
t943 = sin(qJ(2,1));
t1028 = t1027 * t943;
t938 = sin(qJ(3,3));
t1017 = pkin(2) * t938;
t944 = cos(qJ(3,3));
t1014 = pkin(2) * t944;
t945 = cos(qJ(2,3));
t908 = -pkin(5) * t945 + t939 * t1014;
t893 = t934 * t1017 + t908 * t932;
t1021 = 0.1e1 / t893;
t940 = sin(qJ(3,2));
t1016 = pkin(2) * t940;
t946 = cos(qJ(3,2));
t1013 = pkin(2) * t946;
t947 = cos(qJ(2,2));
t909 = -pkin(5) * t947 + t941 * t1013;
t894 = t934 * t1016 + t909 * t932;
t1020 = 0.1e1 / t894;
t942 = sin(qJ(3,1));
t1015 = pkin(2) * t942;
t948 = cos(qJ(3,1));
t1012 = pkin(2) * t948;
t949 = cos(qJ(2,1));
t910 = -pkin(5) * t949 + t943 * t1012;
t895 = t934 * t1015 + t910 * t932;
t1019 = 0.1e1 / t895;
t1018 = m(3) / pkin(2);
t1011 = g(3) * t933;
t1010 = t1021 / t944;
t1009 = t1020 / t946;
t1008 = t1019 / t948;
t1007 = t1021 * t902;
t1006 = t1020 * t903;
t1005 = t1019 * t904;
t1003 = t902 * t934;
t1001 = t903 * t934;
t999 = t904 * t934;
t995 = rSges(3,2) * t920 * t932;
t994 = t931 * t945;
t993 = t931 * t947;
t992 = t931 * t949;
t991 = t932 * t933;
t990 = t934 * t939;
t989 = t934 * t941;
t988 = t934 * t943;
t987 = t934 * t945;
t986 = t934 * t947;
t985 = t934 * t949;
t984 = t938 * t945;
t983 = t940 * t947;
t982 = t942 * t949;
t969 = t905 * t931 + t1011;
t954 = t969 * t945 + t1030;
t981 = (((t968 * t932 - t1003) * rSges(3,1) + t954 * rSges(3,2)) * t944 + t938 * (t995 + (-t905 * t991 + t1003) * rSges(3,2) + t954 * rSges(3,1))) * t1010;
t967 = t906 * t931 + t1011;
t953 = t967 * t947 + t1029;
t980 = (((t966 * t932 - t1001) * rSges(3,1) + t953 * rSges(3,2)) * t946 + t940 * (t995 + (-t906 * t991 + t1001) * rSges(3,2) + t953 * rSges(3,1))) * t1009;
t965 = t907 * t931 + t1011;
t952 = t965 * t949 + t1028;
t979 = (((t964 * t932 - t999) * rSges(3,1) + t952 * rSges(3,2)) * t948 + t942 * (t995 + (-t907 * t991 + t999) * rSges(3,2) + t952 * rSges(3,1))) * t1008;
t918 = m(2) * rSges(2,2) - m(3) * rSges(3,3);
t914 = t918 * t1011;
t950 = m(2) * rSges(2,1);
t978 = (t914 * t945 + (t905 * t994 + t1030) * t918 + (-t1025 * t945 + t939 * t969) * (t950 + (rSges(3,1) * t944 - rSges(3,2) * t938) * m(3))) * t1010;
t977 = (t914 * t947 + (t906 * t993 + t1029) * t918 + (-t1026 * t947 + t941 * t967) * (t950 + (rSges(3,1) * t946 - rSges(3,2) * t940) * m(3))) * t1009;
t976 = (t914 * t949 + (t907 * t992 + t1028) * t918 + (-t1027 * t949 + t943 * t965) * (t950 + (rSges(3,1) * t948 - rSges(3,2) * t942) * m(3))) * t1008;
t975 = ((t931 * t987 + t933 * t939) * t1014 + (t931 * t990 - t933 * t945) * pkin(5)) * t981;
t974 = ((t931 * t986 + t933 * t941) * t1013 + (t931 * t989 - t933 * t947) * pkin(5)) * t980;
t973 = ((t931 * t985 + t933 * t943) * t1012 + (t931 * t988 - t933 * t949) * pkin(5)) * t979;
t957 = t932 * t944 + t938 * t990;
t972 = (t957 * t931 - t933 * t984) * t978;
t956 = t932 * t946 + t940 * t989;
t971 = (t956 * t931 - t933 * t983) * t977;
t955 = t932 * t948 + t942 * t988;
t970 = (t955 * t931 - t933 * t982) * t976;
t960 = t932 * t1017 - t908 * t934;
t959 = t932 * t1016 - t909 * t934;
t958 = t932 * t1015 - t910 * t934;
t927 = m(1) + m(2) + m(3);
t913 = pkin(5) * t943 + t949 * t1012;
t912 = pkin(5) * t941 + t947 * t1013;
t911 = pkin(5) * t939 + t945 * t1014;
t880 = -t913 * t931 + t958 * t933;
t879 = -t912 * t931 + t959 * t933;
t878 = -t911 * t931 + t960 * t933;
t1 = [-t924 * t972 - t925 * t971 - t926 * t970 - m(4) * g(1) + (-(-t880 * t926 + t895 * t923) * t1005 - (-t879 * t925 + t894 * t922) * t1006 - (-t878 * t924 + t893 * t921) * t1007) * t927 + (-t924 * t975 - t925 * t974 - t926 * t973) * t1018; t921 * t972 + t922 * t971 + t923 * t970 - m(4) * g(2) + (-(t880 * t923 + t895 * t926) * t1005 - (t879 * t922 + t894 * t925) * t1006 - (t878 * t921 + t893 * t924) * t1007) * t927 + (t921 * t975 + t922 * t974 + t923 * t973) * t1018; (-t931 * t982 - t955 * t933) * t976 + (-t931 * t983 - t956 * t933) * t977 + (-t931 * t984 - t957 * t933) * t978 - m(4) * g(3) + (-(t913 * t933 + t958 * t931) * t1005 - (t912 * t933 + t959 * t931) * t1006 - (t911 * t933 + t960 * t931) * t1007) * t927 + ((-(-t931 * t943 + t933 * t985) * t1012 - pkin(5) * (t933 * t988 + t992)) * t979 + (-(-t931 * t941 + t933 * t986) * t1013 - pkin(5) * (t933 * t989 + t993)) * t980 + (-(-t931 * t939 + t933 * t987) * t1014 - pkin(5) * (t933 * t990 + t994)) * t981) * t1018;];
taugX  = t1;
