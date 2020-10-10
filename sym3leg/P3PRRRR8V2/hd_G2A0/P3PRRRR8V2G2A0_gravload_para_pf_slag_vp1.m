% Calculate Gravitation load for parallel robot
% P3PRRRR8V2G2A0
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
% Datum: 2020-08-06 17:50
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3PRRRR8V2G2A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(8,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V2G2A0_gravload_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V2G2A0_gravload_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P3PRRRR8V2G2A0_gravload_para_pf_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V2G2A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR8V2G2A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRRR8V2G2A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V2G2A0_gravload_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V2G2A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 17:49:15
% EndTime: 2020-08-06 17:49:16
% DurationCPUTime: 1.38s
% Computational Cost: add. (927->176), mult. (1986->336), div. (36->4), fcn. (1527->22), ass. (0->132)
t984 = legFrame(1,2);
t970 = sin(t984);
t973 = cos(t984);
t957 = t973 * g(1) - t970 * g(2);
t978 = sin(pkin(8));
t967 = g(3) * t978;
t980 = cos(pkin(8));
t1080 = -t957 * t980 + t967;
t954 = t970 * g(1) + t973 * g(2);
t979 = sin(pkin(4));
t981 = cos(pkin(4));
t1007 = t1080 * t981 - t954 * t979;
t990 = sin(qJ(2,1));
t1083 = t1007 * t990;
t983 = legFrame(2,2);
t969 = sin(t983);
t972 = cos(t983);
t956 = t972 * g(1) - t969 * g(2);
t1079 = -t956 * t980 + t967;
t953 = t969 * g(1) + t972 * g(2);
t1008 = t1079 * t981 - t953 * t979;
t988 = sin(qJ(2,2));
t1082 = t1008 * t988;
t982 = legFrame(3,2);
t968 = sin(t982);
t971 = cos(t982);
t955 = t971 * g(1) - t968 * g(2);
t1078 = -t955 * t980 + t967;
t952 = t968 * g(1) + t971 * g(2);
t1009 = t1078 * t981 - t952 * t979;
t986 = sin(qJ(2,3));
t1081 = t1009 * t986;
t985 = sin(qJ(3,3));
t1077 = pkin(2) * t985;
t987 = sin(qJ(3,2));
t1076 = pkin(2) * t987;
t989 = sin(qJ(3,1));
t1075 = pkin(2) * t989;
t991 = cos(qJ(3,3));
t1074 = pkin(3) * t991 ^ 2;
t993 = cos(qJ(3,2));
t1073 = pkin(3) * t993 ^ 2;
t995 = cos(qJ(3,1));
t1072 = pkin(3) * t995 ^ 2;
t1071 = pkin(3) * t991;
t1070 = pkin(3) * t993;
t1069 = pkin(3) * t995;
t1068 = g(3) * t980;
t1067 = m(3) / pkin(3);
t1015 = t955 * t978 + t1068;
t992 = cos(qJ(2,3));
t1003 = t1015 * t992 - t1081;
t1042 = t979 * t980;
t1048 = rSges(3,2) * t967 * t979;
t1056 = t952 * t981;
t1031 = t981 * t985;
t1040 = t979 * t986;
t998 = pkin(7) + pkin(6);
t959 = pkin(2) * t986 - t998 * t992;
t937 = pkin(3) * t1031 + t959 * t979;
t928 = 0.1e1 / (pkin(2) * t1031 + t1040 * t1074 + t937 * t991);
t1066 = (((-t1078 * t979 - t1056) * rSges(3,1) + t1003 * rSges(3,2)) * t991 + t985 * (t1048 + (-t955 * t1042 + t1056) * rSges(3,2) + t1003 * rSges(3,1))) * t928;
t1013 = t956 * t978 + t1068;
t994 = cos(qJ(2,2));
t1002 = t1013 * t994 - t1082;
t1054 = t953 * t981;
t1029 = t981 * t987;
t1038 = t979 * t988;
t960 = pkin(2) * t988 - t998 * t994;
t938 = pkin(3) * t1029 + t960 * t979;
t929 = 0.1e1 / (pkin(2) * t1029 + t1038 * t1073 + t938 * t993);
t1065 = (((-t1079 * t979 - t1054) * rSges(3,1) + t1002 * rSges(3,2)) * t993 + t987 * (t1048 + (-t956 * t1042 + t1054) * rSges(3,2) + t1002 * rSges(3,1))) * t929;
t1011 = t957 * t978 + t1068;
t996 = cos(qJ(2,1));
t1001 = t1011 * t996 - t1083;
t1052 = t954 * t981;
t1027 = t981 * t989;
t1036 = t979 * t990;
t961 = pkin(2) * t990 - t998 * t996;
t939 = pkin(3) * t1027 + t961 * t979;
t930 = 0.1e1 / (pkin(2) * t1027 + t1036 * t1072 + t939 * t995);
t1064 = (((-t1080 * t979 - t1052) * rSges(3,1) + t1001 * rSges(3,2)) * t995 + t989 * (t1048 + (-t957 * t1042 + t1052) * rSges(3,2) + t1001 * rSges(3,1))) * t930;
t1022 = m(2) * rSges(2,1) + pkin(2) * m(3);
t1045 = t978 * t992;
t965 = (-pkin(6) - rSges(3,3)) * m(3) + m(2) * rSges(2,2);
t958 = t965 * t1068;
t1063 = (t958 * t992 + (t955 * t1045 - t1081) * t965 + (t1009 * t992 + t1015 * t986) * ((rSges(3,1) * t991 - rSges(3,2) * t985) * m(3) + t1022)) * t928;
t1044 = t978 * t994;
t1062 = (t958 * t994 + (t956 * t1044 - t1082) * t965 + (t1008 * t994 + t1013 * t988) * ((rSges(3,1) * t993 - rSges(3,2) * t987) * m(3) + t1022)) * t929;
t1043 = t978 * t996;
t1061 = (t958 * t996 + (t957 * t1043 - t1083) * t965 + (t1007 * t996 + t1011 * t990) * ((rSges(3,1) * t995 - rSges(3,2) * t989) * m(3) + t1022)) * t930;
t1060 = t928 * t952;
t1059 = t929 * t953;
t1058 = t930 * t954;
t1047 = t978 * t979;
t1046 = t978 * t981;
t1041 = t979 * t985;
t1039 = t979 * t987;
t1037 = t979 * t989;
t1035 = t980 * t981;
t1034 = t980 * t991;
t1033 = t980 * t993;
t1032 = t980 * t995;
t1030 = t981 * t986;
t1028 = t981 * t988;
t1026 = t981 * t990;
t1025 = t981 * t992;
t1024 = t981 * t994;
t1023 = t981 * t996;
t962 = pkin(2) * t992 + t986 * t998;
t1021 = ((t978 * t1025 + t980 * t986) * t1071 + t962 * t1046 + t959 * t980) * t1066;
t963 = pkin(2) * t994 + t988 * t998;
t1020 = ((t978 * t1024 + t980 * t988) * t1070 + t963 * t1046 + t960 * t980) * t1065;
t964 = pkin(2) * t996 + t990 * t998;
t1019 = ((t978 * t1023 + t980 * t990) * t1069 + t964 * t1046 + t961 * t980) * t1064;
t946 = t978 * t1030 - t980 * t992;
t1018 = (t991 * t1047 + t985 * t946) * t1063;
t947 = t978 * t1028 - t980 * t994;
t1017 = (t993 * t1047 + t987 * t947) * t1062;
t948 = t978 * t1026 - t980 * t996;
t1016 = (t995 * t1047 + t989 * t948) * t1061;
t1006 = pkin(3) * t1041 - t959 * t981;
t1005 = pkin(3) * t1039 - t960 * t981;
t1004 = pkin(3) * t1037 - t961 * t981;
t974 = m(1) + m(2) + m(3);
t951 = t980 * t1026 + t1043;
t950 = t980 * t1028 + t1044;
t949 = t980 * t1030 + t1045;
t933 = t1004 * t980 - t978 * t964;
t932 = t1005 * t980 - t978 * t963;
t931 = t1006 * t980 - t978 * t962;
t1 = [-t971 * t1018 - t972 * t1017 - t973 * t1016 - m(4) * g(1) + (-((t970 * t1036 + t951 * t973) * t1072 + (-t933 * t973 + t939 * t970) * t995 + (-t973 * t1042 + t981 * t970) * t1075) * t1058 - ((t969 * t1038 + t950 * t972) * t1073 + (-t932 * t972 + t938 * t969) * t993 + (-t972 * t1042 + t981 * t969) * t1076) * t1059 - ((t968 * t1040 + t949 * t971) * t1074 + (-t931 * t971 + t937 * t968) * t991 + (-t971 * t1042 + t981 * t968) * t1077) * t1060) * t974 + (-t973 * t1019 - t972 * t1020 - t971 * t1021) * t1067; t968 * t1018 + t969 * t1017 + t970 * t1016 - m(4) * g(2) + (-(-(-t973 * t1036 + t951 * t970) * t1072 + (t933 * t970 + t973 * t939) * t995 + (t970 * t1042 + t973 * t981) * t1075) * t1058 - (-(-t972 * t1038 + t950 * t969) * t1073 + (t932 * t969 + t972 * t938) * t993 + (t969 * t1042 + t972 * t981) * t1076) * t1059 - (-(-t971 * t1040 + t949 * t968) * t1074 + (t931 * t968 + t971 * t937) * t991 + (t968 * t1042 + t971 * t981) * t1077) * t1060) * t974 + (t970 * t1019 + t969 * t1020 + t968 * t1021) * t1067; (-t979 * t1032 - t989 * t951) * t1061 + (-t979 * t1033 - t987 * t950) * t1062 + (-t979 * t1034 - t985 * t949) * t1063 - m(4) * g(3) + (-(-t948 * t1072 + t964 * t1032 + (pkin(2) * t1037 + t1004 * t995) * t978) * t1058 - (-t947 * t1073 + t963 * t1033 + (pkin(2) * t1039 + t1005 * t993) * t978) * t1059 - (-t946 * t1074 + t962 * t1034 + (pkin(2) * t1041 + t1006 * t991) * t978) * t1060) * t974 + (((-t980 * t1023 + t978 * t990) * t1069 - t964 * t1035 + t978 * t961) * t1064 + ((-t980 * t1024 + t978 * t988) * t1070 - t963 * t1035 + t978 * t960) * t1065 + ((-t980 * t1025 + t978 * t986) * t1071 - t962 * t1035 + t978 * t959) * t1066) * t1067;];
taugX  = t1;
