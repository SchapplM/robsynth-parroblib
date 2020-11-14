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
% Datum: 2020-08-06 17:36
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3PRRRR8V2G1A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(8,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V2G1A0_gravload_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V2G1A0_gravload_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P3PRRRR8V2G1A0_gravload_para_pf_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V2G1A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR8V2G1A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRRR8V2G1A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V2G1A0_gravload_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V2G1A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 17:35:20
% EndTime: 2020-08-06 17:35:22
% DurationCPUTime: 1.38s
% Computational Cost: add. (705->159), mult. (1439->303), div. (24->7), fcn. (1170->22), ass. (0->127)
t1004 = cos(qJ(2,1));
t992 = legFrame(1,3);
t978 = sin(t992);
t981 = cos(t992);
t959 = -t978 * g(1) + t981 * g(2);
t962 = t981 * g(1) + t978 * g(2);
t986 = sin(pkin(8));
t988 = cos(pkin(8));
t1019 = t959 * t986 + t962 * t988;
t1020 = t959 * t988 - t962 * t986;
t987 = sin(pkin(4));
t1063 = g(3) * t987;
t989 = cos(pkin(4));
t1071 = t1020 * t989 + t1063;
t998 = sin(qJ(2,1));
t1010 = t1019 * t1004 + t1071 * t998;
t1002 = cos(qJ(2,2));
t991 = legFrame(2,3);
t977 = sin(t991);
t980 = cos(t991);
t958 = -t977 * g(1) + t980 * g(2);
t961 = t980 * g(1) + t977 * g(2);
t1021 = t958 * t986 + t961 * t988;
t1022 = t958 * t988 - t961 * t986;
t1070 = t1022 * t989 + t1063;
t996 = sin(qJ(2,2));
t1011 = t1021 * t1002 + t1070 * t996;
t1000 = cos(qJ(2,3));
t990 = legFrame(3,3);
t976 = sin(t990);
t979 = cos(t990);
t957 = -t976 * g(1) + t979 * g(2);
t960 = t979 * g(1) + t976 * g(2);
t1023 = t957 * t986 + t960 * t988;
t1024 = t957 * t988 - t960 * t986;
t1069 = t1024 * t989 + t1063;
t994 = sin(qJ(2,3));
t1012 = t1023 * t1000 + t1069 * t994;
t1062 = g(3) * t989;
t997 = sin(qJ(3,1));
t1045 = t987 * t997;
t1007 = pkin(7) + pkin(6);
t975 = t1007 * t1004;
t965 = pkin(2) * t998 - t975;
t1016 = pkin(3) * t1045 - t965 * t989;
t1032 = t998 * t1007;
t968 = pkin(2) * t1004 + t1032;
t1074 = t1016 * t986 + t968 * t988;
t995 = sin(qJ(3,2));
t1046 = t987 * t995;
t974 = t1007 * t1002;
t964 = pkin(2) * t996 - t974;
t1017 = pkin(3) * t1046 - t964 * t989;
t1033 = t996 * t1007;
t967 = pkin(2) * t1002 + t1033;
t1073 = t1017 * t986 + t967 * t988;
t993 = sin(qJ(3,3));
t1047 = t987 * t993;
t973 = t1007 * t1000;
t963 = pkin(2) * t994 - t973;
t1018 = pkin(3) * t1047 - t963 * t989;
t1034 = t994 * t1007;
t966 = pkin(2) * t1000 + t1034;
t1072 = t1018 * t986 + t966 * t988;
t999 = cos(qJ(3,3));
t983 = t999 ^ 2;
t1068 = pkin(3) * t983;
t1001 = cos(qJ(3,2));
t984 = t1001 ^ 2;
t1067 = pkin(3) * t984;
t1003 = cos(qJ(3,1));
t985 = t1003 ^ 2;
t1066 = pkin(3) * t985;
t1065 = pkin(3) * t987;
t982 = m(1) + m(2) + m(3);
t1064 = g(3) * t982;
t1061 = m(3) / pkin(3);
t1060 = rSges(3,2) * t987;
t1037 = rSges(3,2) * t1062;
t1043 = t989 * t993;
t1044 = t987 * t999;
t970 = t999 * pkin(3) + pkin(2);
t954 = t994 * t970 - t973;
t1059 = (((t1024 * t987 - t1062) * rSges(3,1) + t1012 * rSges(3,2)) * t999 + t993 * (t1012 * rSges(3,1) - t1024 * t1060 + t1037)) / (t970 * t1043 + t954 * t1044);
t1036 = t987 * t1001;
t1041 = t989 * t995;
t971 = t1001 * pkin(3) + pkin(2);
t955 = t996 * t971 - t974;
t1058 = (((t1022 * t987 - t1062) * rSges(3,1) + t1011 * rSges(3,2)) * t1001 + t995 * (t1011 * rSges(3,1) - t1022 * t1060 + t1037)) / (t955 * t1036 + t971 * t1041);
t1035 = t987 * t1003;
t1039 = t989 * t997;
t972 = t1003 * pkin(3) + pkin(2);
t956 = t998 * t972 - t975;
t1057 = (((t1020 * t987 - t1062) * rSges(3,1) + t1010 * rSges(3,2)) * t1003 + t997 * (t1010 * rSges(3,1) - t1020 * t1060 + t1037)) / (t956 * t1035 + t972 * t1039);
t1031 = m(2) * rSges(2,1) + pkin(2) * m(3);
t924 = 0.1e1 / (t994 * t983 * t1065 + (pkin(3) * t1043 + t963 * t987) * t999 + pkin(2) * t1043);
t969 = (-pkin(6) - rSges(3,3)) * m(3) + m(2) * rSges(2,2);
t1056 = (t1012 * t969 + (-t1000 * t1069 + t1023 * t994) * ((rSges(3,1) * t999 - rSges(3,2) * t993) * m(3) + t1031)) * t924;
t925 = 0.1e1 / (t996 * t984 * t1065 + (pkin(3) * t1041 + t964 * t987) * t1001 + pkin(2) * t1041);
t1055 = (t1011 * t969 + (-t1002 * t1070 + t1021 * t996) * ((rSges(3,1) * t1001 - rSges(3,2) * t995) * m(3) + t1031)) * t925;
t926 = 0.1e1 / (t998 * t985 * t1065 + (pkin(3) * t1039 + t965 * t987) * t1003 + pkin(2) * t1039);
t1054 = (t1010 * t969 + (-t1004 * t1071 + t1019 * t998) * ((rSges(3,1) * t1003 - rSges(3,2) * t997) * m(3) + t1031)) * t926;
t1053 = (t970 * t1000 + t1034) * t989;
t1052 = (t971 * t1002 + t1033) * t989;
t1051 = (t972 * t1004 + t1032) * t989;
t1042 = t989 * t994;
t1040 = t989 * t996;
t1038 = t989 * t998;
t1030 = pkin(2) * t1047;
t1029 = pkin(2) * t1046;
t1028 = pkin(2) * t1045;
t950 = t986 * t1004 + t988 * t1038;
t949 = t986 * t1002 + t988 * t1040;
t948 = t986 * t1000 + t988 * t1042;
t947 = -t988 * t1004 + t986 * t1038;
t946 = -t988 * t1002 + t986 * t1040;
t945 = -t988 * t1000 + t986 * t1042;
t941 = t988 * t978 + t981 * t986;
t940 = t988 * t977 + t980 * t986;
t939 = t988 * t976 + t979 * t986;
t938 = -t986 * t978 + t981 * t988;
t937 = -t986 * t977 + t980 * t988;
t936 = -t986 * t976 + t979 * t988;
t929 = -t1016 * t988 + t986 * t968;
t928 = -t1017 * t988 + t986 * t967;
t927 = -t1018 * t988 + t986 * t966;
t1 = [(-t938 * t1035 - (t1004 * t941 + t938 * t1038) * t997) * t1054 + (-t937 * t1036 - (t1002 * t940 + t937 * t1040) * t995) * t1055 + (-t936 * t1044 - (t1000 * t939 + t936 * t1042) * t993) * t1056 - m(4) * g(1) + (-(-(t947 * t981 + t978 * t950) * t1066 + (t1074 * t981 - t929 * t978) * t1003 + t941 * t1028) * t926 - (-(t946 * t980 + t977 * t949) * t1067 + (t1073 * t980 - t928 * t977) * t1001 + t940 * t1029) * t925 - (-(t945 * t979 + t976 * t948) * t1068 + (t1072 * t979 - t927 * t976) * t999 + t939 * t1030) * t924) * t1064 + ((-t938 * t1051 + t956 * t941) * t1057 + (-t937 * t1052 + t955 * t940) * t1058 + (-t936 * t1053 + t954 * t939) * t1059) * t1061; (-t941 * t1035 - (-t1004 * t938 + t941 * t1038) * t997) * t1054 + (-t940 * t1036 - (-t1002 * t937 + t940 * t1040) * t995) * t1055 + (-t939 * t1044 - (-t1000 * t936 + t939 * t1042) * t993) * t1056 - m(4) * g(2) + (-((-t978 * t947 + t950 * t981) * t1066 + (t1074 * t978 + t929 * t981) * t1003 - t938 * t1028) * t926 - ((-t977 * t946 + t949 * t980) * t1067 + (t1073 * t977 + t928 * t980) * t1001 - t937 * t1029) * t925 - ((-t976 * t945 + t948 * t979) * t1068 + (t1072 * t976 + t927 * t979) * t999 - t936 * t1030) * t924) * t1064 + ((-t941 * t1051 - t956 * t938) * t1057 + (-t940 * t1052 - t955 * t937) * t1058 + (-t939 * t1053 - t954 * t936) * t1059) * t1061; (-m(4) - 0.3e1 * t982) * g(3);];
taugX  = t1;
