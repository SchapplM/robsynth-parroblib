% Calculate Gravitation load for parallel robot
% P3PRRRR8V2G1A0
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [3x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 17:36
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3PRRRR8V2G1A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(8,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V2G1A0_gravload_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V2G1A0_gravload_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P3PRRRR8V2G1A0_gravload_para_pf_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V2G1A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR8V2G1A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRRR8V2G1A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V2G1A0_gravload_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V2G1A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 17:35:36
% EndTime: 2020-08-06 17:35:37
% DurationCPUTime: 1.30s
% Computational Cost: add. (705->159), mult. (1337->298), div. (24->7), fcn. (1170->22), ass. (0->127)
t971 = legFrame(3,3);
t957 = sin(t971);
t960 = cos(t971);
t938 = -t957 * g(1) + t960 * g(2);
t941 = t960 * g(1) + t957 * g(2);
t967 = sin(pkin(8));
t969 = cos(pkin(8));
t1000 = t938 * t967 + t941 * t969;
t1001 = t938 * t969 - t941 * t967;
t968 = sin(pkin(4));
t1039 = g(3) * t968;
t970 = cos(pkin(4));
t1045 = t1001 * t970 + t1039;
t975 = sin(qJ(2,3));
t981 = cos(qJ(2,3));
t992 = t1000 * t981 + t1045 * t975;
t972 = legFrame(2,3);
t958 = sin(t972);
t961 = cos(t972);
t939 = -t958 * g(1) + t961 * g(2);
t942 = t961 * g(1) + t958 * g(2);
t999 = t939 * t969 - t942 * t967;
t1046 = t999 * t970 + t1039;
t977 = sin(qJ(2,2));
t983 = cos(qJ(2,2));
t998 = t939 * t967 + t942 * t969;
t991 = t1046 * t977 + t998 * t983;
t973 = legFrame(1,3);
t959 = sin(t973);
t962 = cos(t973);
t940 = -t959 * g(1) + t962 * g(2);
t943 = t962 * g(1) + t959 * g(2);
t997 = t940 * t969 - t943 * t967;
t1047 = t997 * t970 + t1039;
t979 = sin(qJ(2,1));
t985 = cos(qJ(2,1));
t996 = t940 * t967 + t943 * t969;
t990 = t1047 * t979 + t996 * t985;
t1038 = g(3) * t970;
t987 = pkin(7) + pkin(6);
t1009 = t979 * t987;
t949 = pkin(2) * t985 + t1009;
t978 = sin(qJ(3,1));
t1021 = t968 * t978;
t955 = t987 * t985;
t946 = pkin(2) * t979 - t955;
t993 = pkin(3) * t1021 - t946 * t970;
t1050 = t949 * t969 + t993 * t967;
t1010 = t977 * t987;
t948 = pkin(2) * t983 + t1010;
t976 = sin(qJ(3,2));
t1022 = t968 * t976;
t954 = t987 * t983;
t945 = pkin(2) * t977 - t954;
t994 = pkin(3) * t1022 - t945 * t970;
t1049 = t948 * t969 + t994 * t967;
t1011 = t975 * t987;
t947 = pkin(2) * t981 + t1011;
t974 = sin(qJ(3,3));
t1023 = t968 * t974;
t953 = t987 * t981;
t944 = pkin(2) * t975 - t953;
t995 = pkin(3) * t1023 - t944 * t970;
t1048 = t947 * t969 + t995 * t967;
t980 = cos(qJ(3,3));
t964 = t980 ^ 2;
t1044 = pkin(3) * t964;
t982 = cos(qJ(3,2));
t965 = t982 ^ 2;
t1043 = pkin(3) * t965;
t984 = cos(qJ(3,1));
t966 = t984 ^ 2;
t1042 = pkin(3) * t966;
t1041 = pkin(3) * t968;
t963 = m(1) + m(2) + m(3);
t1040 = g(3) * t963;
t1037 = mrSges(3,2) * t968;
t1036 = m(3) * pkin(2) + mrSges(2,1);
t1008 = mrSges(3,2) * t1038;
t1017 = t970 * t974;
t1020 = t968 * t980;
t950 = t980 * pkin(3) + pkin(2);
t932 = t975 * t950 - t953;
t1035 = (((t1001 * t968 - t1038) * mrSges(3,1) + t992 * mrSges(3,2)) * t980 + (t992 * mrSges(3,1) - t1001 * t1037 + t1008) * t974) / (t950 * t1017 + t932 * t1020);
t1015 = t970 * t976;
t1019 = t968 * t982;
t951 = t982 * pkin(3) + pkin(2);
t933 = t977 * t951 - t954;
t1034 = (((t999 * t968 - t1038) * mrSges(3,1) + t991 * mrSges(3,2)) * t982 + (t991 * mrSges(3,1) - t999 * t1037 + t1008) * t976) / (t951 * t1015 + t933 * t1019);
t1013 = t970 * t978;
t1018 = t968 * t984;
t952 = t984 * pkin(3) + pkin(2);
t934 = t979 * t952 - t955;
t1033 = (((t997 * t968 - t1038) * mrSges(3,1) + t990 * mrSges(3,2)) * t984 + (t990 * mrSges(3,1) - t997 * t1037 + t1008) * t978) / (t952 * t1013 + t934 * t1018);
t905 = 0.1e1 / (t975 * t964 * t1041 + (pkin(3) * t1017 + t944 * t968) * t980 + pkin(2) * t1017);
t956 = m(3) * pkin(6) - mrSges(2,2) + mrSges(3,3);
t1032 = (-t992 * t956 + (t1000 * t975 - t1045 * t981) * (t980 * mrSges(3,1) - mrSges(3,2) * t974 + t1036)) * t905;
t906 = 0.1e1 / (t977 * t965 * t1041 + (pkin(3) * t1015 + t945 * t968) * t982 + pkin(2) * t1015);
t1031 = (-t991 * t956 + (-t1046 * t983 + t998 * t977) * (t982 * mrSges(3,1) - mrSges(3,2) * t976 + t1036)) * t906;
t907 = 0.1e1 / (t979 * t966 * t1041 + (pkin(3) * t1013 + t946 * t968) * t984 + pkin(2) * t1013);
t1030 = (-t990 * t956 + (-t1047 * t985 + t996 * t979) * (t984 * mrSges(3,1) - mrSges(3,2) * t978 + t1036)) * t907;
t1029 = (t950 * t981 + t1011) * t970;
t1028 = (t951 * t983 + t1010) * t970;
t1027 = (t952 * t985 + t1009) * t970;
t1016 = t970 * t975;
t1014 = t970 * t977;
t1012 = t970 * t979;
t1007 = pkin(2) * t1023;
t1006 = pkin(2) * t1022;
t1005 = pkin(2) * t1021;
t989 = 0.1e1 / pkin(3);
t928 = t969 * t1012 + t967 * t985;
t927 = t969 * t1014 + t967 * t983;
t926 = t969 * t1016 + t967 * t981;
t925 = t967 * t1012 - t969 * t985;
t924 = t967 * t1014 - t969 * t983;
t923 = t967 * t1016 - t969 * t981;
t922 = t969 * t959 + t962 * t967;
t921 = t969 * t958 + t961 * t967;
t920 = t969 * t957 + t960 * t967;
t919 = -t967 * t959 + t962 * t969;
t918 = -t967 * t958 + t961 * t969;
t917 = -t967 * t957 + t960 * t969;
t910 = t967 * t949 - t993 * t969;
t909 = t967 * t948 - t994 * t969;
t908 = t967 * t947 - t995 * t969;
t1 = [(-t919 * t1018 - (t919 * t1012 + t985 * t922) * t978) * t1030 + (-t918 * t1019 - (t918 * t1014 + t983 * t921) * t976) * t1031 + (-t917 * t1020 - (t917 * t1016 + t981 * t920) * t974) * t1032 - g(1) * m(4) + ((-t919 * t1027 + t934 * t922) * t1033 + (-t918 * t1028 + t933 * t921) * t1034 + (-t917 * t1029 + t932 * t920) * t1035) * t989 + (-(-(t925 * t962 + t959 * t928) * t1042 + (t1050 * t962 - t910 * t959) * t984 + t922 * t1005) * t907 - (-(t924 * t961 + t958 * t927) * t1043 + (t1049 * t961 - t909 * t958) * t982 + t921 * t1006) * t906 - (-(t923 * t960 + t957 * t926) * t1044 + (t1048 * t960 - t908 * t957) * t980 + t920 * t1007) * t905) * t1040; (-t922 * t1018 - (t922 * t1012 - t985 * t919) * t978) * t1030 + (-t921 * t1019 - (t921 * t1014 - t983 * t918) * t976) * t1031 + (-t920 * t1020 - (t920 * t1016 - t981 * t917) * t974) * t1032 - g(2) * m(4) + ((-t922 * t1027 - t934 * t919) * t1033 + (-t921 * t1028 - t933 * t918) * t1034 + (-t920 * t1029 - t932 * t917) * t1035) * t989 + (-((-t959 * t925 + t928 * t962) * t1042 + (t1050 * t959 + t910 * t962) * t984 - t919 * t1005) * t907 - ((-t958 * t924 + t927 * t961) * t1043 + (t1049 * t958 + t909 * t961) * t982 - t918 * t1006) * t906 - ((-t957 * t923 + t926 * t960) * t1044 + (t1048 * t957 + t908 * t960) * t980 - t917 * t1007) * t905) * t1040; (-m(4) - 0.3e1 * t963) * g(3);];
taugX  = t1;
