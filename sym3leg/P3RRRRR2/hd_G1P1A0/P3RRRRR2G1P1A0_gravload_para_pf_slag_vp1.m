% Calculate Gravitation load for parallel robot
% P3RRRRR2G1P1A0
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% Datum: 2020-03-09 21:05
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRRRR2G1P1A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR2G1P1A0_gravload_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR2G1P1A0_gravload_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3RRRRR2G1P1A0_gravload_para_pf_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRRRR2G1P1A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRRRR2G1P1A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRRRR2G1P1A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR2G1P1A0_gravload_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR2G1P1A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:04:10
% EndTime: 2020-03-09 21:04:11
% DurationCPUTime: 0.75s
% Computational Cost: add. (783->152), mult. (1383->210), div. (57->14), fcn. (771->72), ass. (0->111)
t1041 = sin(qJ(3,1));
t1044 = cos(qJ(3,1));
t1092 = t1044 * rSges(3,1) - t1041 * rSges(3,2);
t1040 = sin(qJ(3,2));
t1043 = cos(qJ(3,2));
t1091 = t1043 * rSges(3,1) - t1040 * rSges(3,2);
t1039 = sin(qJ(3,3));
t1042 = cos(qJ(3,3));
t1090 = t1042 * rSges(3,1) - t1039 * rSges(3,2);
t1089 = -2 * pkin(1);
t1088 = -m(3) / 0.2e1;
t1087 = m(3) / 0.2e1;
t1086 = m(1) * rSges(1,2);
t1036 = legFrame(3,3);
t1015 = sin(t1036);
t1018 = cos(t1036);
t984 = -t1015 * g(1) + t1018 * g(2);
t1084 = rSges(3,1) * t984;
t1037 = legFrame(2,3);
t1016 = sin(t1037);
t1019 = cos(t1037);
t985 = -t1016 * g(1) + t1019 * g(2);
t1083 = rSges(3,1) * t985;
t1038 = legFrame(1,3);
t1017 = sin(t1038);
t1020 = cos(t1038);
t986 = -t1017 * g(1) + t1020 * g(2);
t1082 = rSges(3,1) * t986;
t987 = t1018 * g(1) + t1015 * g(2);
t1081 = rSges(3,1) * t987;
t988 = t1019 * g(1) + t1016 * g(2);
t1080 = rSges(3,1) * t988;
t989 = t1020 * g(1) + t1017 * g(2);
t1079 = rSges(3,1) * t989;
t1078 = rSges(3,3) * t984;
t1077 = rSges(3,3) * t985;
t1076 = rSges(3,3) * t986;
t1075 = rSges(3,3) * t987;
t1074 = rSges(3,3) * t988;
t1073 = rSges(3,3) * t989;
t1066 = qJ(2,1) - qJ(3,1);
t1065 = qJ(2,1) + qJ(3,1);
t1064 = qJ(2,2) - qJ(3,2);
t1063 = qJ(2,2) + qJ(3,2);
t1062 = qJ(2,3) - qJ(3,3);
t1061 = qJ(2,3) + qJ(3,3);
t1008 = qJ(1,1) + t1038;
t1007 = qJ(1,2) + t1037;
t1006 = qJ(1,3) + t1036;
t1046 = 1 / pkin(1);
t1060 = 0.1e1 / sin(qJ(2,3)) * t1046;
t1059 = 0.1e1 / sin(qJ(2,2)) * t1046;
t1058 = 0.1e1 / sin(qJ(2,1)) * t1046;
t1045 = 1 / pkin(2);
t1057 = t1045 * t1046;
t1004 = qJ(2,1) + t1008;
t1003 = qJ(2,2) + t1007;
t1002 = qJ(2,3) + t1006;
t1056 = m(1) * rSges(1,1) + (m(2) + m(3)) * pkin(1);
t1033 = qJ(1,3) + qJ(2,3);
t1009 = sin(t1033);
t1012 = cos(t1033);
t1021 = qJ(1,3) + t1061;
t1022 = qJ(1,3) + t1062;
t972 = m(2) * (-rSges(2,1) * t984 + rSges(2,2) * t987);
t975 = m(2) * (rSges(2,1) * t987 + rSges(2,2) * t984);
t978 = rSges(3,2) * t984;
t981 = rSges(3,2) * t987;
t966 = sin(qJ(1,3)) * (t1056 * t987 + t984 * t1086) + (-m(3) * t1075 + t972) * t1012 + (-m(3) * t1078 + t975) * t1009 + (-t1056 * t984 + t987 * t1086) * cos(qJ(1,3)) + ((t981 + t1084) * cos(t1022) + (t978 - t1081) * sin(t1022)) * t1088 + ((t981 - t1084) * cos(t1021) + (t978 + t1081) * sin(t1021)) * t1087;
t1055 = t966 * t1060;
t1034 = qJ(2,2) + qJ(1,2);
t1010 = sin(t1034);
t1013 = cos(t1034);
t1023 = qJ(1,2) + t1063;
t1024 = qJ(1,2) + t1064;
t973 = m(2) * (-rSges(2,1) * t985 + rSges(2,2) * t988);
t976 = m(2) * (rSges(2,1) * t988 + rSges(2,2) * t985);
t979 = rSges(3,2) * t985;
t982 = rSges(3,2) * t988;
t967 = sin(qJ(1,2)) * (t1056 * t988 + t985 * t1086) + (-m(3) * t1074 + t973) * t1013 + (-m(3) * t1077 + t976) * t1010 + (-t1056 * t985 + t988 * t1086) * cos(qJ(1,2)) + ((t982 + t1083) * cos(t1024) + (t979 - t1080) * sin(t1024)) * t1088 + ((t982 - t1083) * cos(t1023) + (t979 + t1080) * sin(t1023)) * t1087;
t1054 = t967 * t1059;
t1035 = qJ(1,1) + qJ(2,1);
t1011 = sin(t1035);
t1014 = cos(t1035);
t1025 = qJ(1,1) + t1065;
t1026 = qJ(1,1) + t1066;
t974 = m(2) * (-rSges(2,1) * t986 + rSges(2,2) * t989);
t977 = m(2) * (rSges(2,1) * t989 + rSges(2,2) * t986);
t980 = rSges(3,2) * t986;
t983 = rSges(3,2) * t989;
t968 = sin(qJ(1,1)) * (t1056 * t989 + t986 * t1086) + (-m(3) * t1073 + t974) * t1014 + (-m(3) * t1076 + t977) * t1011 + (-t1056 * t986 + t989 * t1086) * cos(qJ(1,1)) + ((t983 + t1082) * cos(t1026) + (t980 - t1079) * sin(t1026)) * t1088 + ((t983 - t1082) * cos(t1025) + (t980 + t1079) * sin(t1025)) * t1087;
t1053 = t968 * t1058;
t1052 = t1039 * t1060;
t1051 = t1040 * t1059;
t1050 = t1041 * t1058;
t969 = t975 * t1009 + t972 * t1012 + ((-t1090 * t984 - t1075) * t1012 + (t1090 * t987 - t1078) * t1009) * m(3);
t1049 = t969 / (sin(t1061) + sin(t1062)) * t1057;
t970 = t976 * t1010 + t973 * t1013 + ((-t1091 * t985 - t1074) * t1013 + (t1091 * t988 - t1077) * t1010) * m(3);
t1048 = t970 / (sin(t1063) + sin(t1064)) * t1057;
t971 = t977 * t1011 + t974 * t1014 + ((-t1092 * t986 - t1073) * t1014 + (t1092 * t989 - t1076) * t1011) * m(3);
t1047 = t971 / (sin(t1065) + sin(t1066)) * t1057;
t1032 = 0.1e1 / t1044;
t1031 = 0.1e1 / t1043;
t1030 = 0.1e1 / t1042;
t1001 = -qJ(3,1) + t1004;
t1000 = qJ(3,1) + t1004;
t999 = -qJ(3,2) + t1003;
t998 = qJ(3,2) + t1003;
t997 = -qJ(3,3) + t1002;
t996 = qJ(3,3) + t1002;
t1 = [cos(t1004) * t1053 + (cos(t1008) * t1089 + (-cos(t1001) - cos(t1000)) * pkin(2)) * t1047 + cos(t1003) * t1054 + (cos(t1007) * t1089 + (-cos(t999) - cos(t998)) * pkin(2)) * t1048 + cos(t1002) * t1055 + (cos(t1006) * t1089 + (-cos(t997) - cos(t996)) * pkin(2)) * t1049 - m(4) * g(1); sin(t1004) * t1053 + (sin(t1008) * t1089 + (-sin(t1001) - sin(t1000)) * pkin(2)) * t1047 + sin(t1003) * t1054 + (sin(t1007) * t1089 + (-sin(t999) - sin(t998)) * pkin(2)) * t1048 + sin(t1002) * t1055 + (sin(t1006) * t1089 + (-sin(t997) - sin(t996)) * pkin(2)) * t1049 - m(4) * g(2); t1030 * t966 * t1052 + t1031 * t967 * t1051 + t1032 * t968 * t1050 - m(4) * g(3) + ((t1032 * (-g(3) * t1092 + (t986 * t1011 + t989 * t1014) * (t1041 * rSges(3,1) + t1044 * rSges(3,2))) + t1031 * (-g(3) * t1091 + (t985 * t1010 + t988 * t1013) * (t1040 * rSges(3,1) + t1043 * rSges(3,2))) + t1030 * (-g(3) * t1090 + (t984 * t1009 + t987 * t1012) * (t1039 * rSges(3,1) + t1042 * rSges(3,2)))) * m(3) - (cos(qJ(2,1)) * pkin(1) + t1044 * pkin(2)) / t1044 ^ 2 * t971 * t1050 - (cos(qJ(2,2)) * pkin(1) + t1043 * pkin(2)) / t1043 ^ 2 * t970 * t1051 - (cos(qJ(2,3)) * pkin(1) + t1042 * pkin(2)) / t1042 ^ 2 * t969 * t1052) * t1045;];
taugX  = t1;
