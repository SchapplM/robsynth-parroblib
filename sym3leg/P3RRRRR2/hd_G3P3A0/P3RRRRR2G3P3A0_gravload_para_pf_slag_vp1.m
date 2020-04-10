% Calculate Gravitation load for parallel robot
% P3RRRRR2G3P3A0
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
% Datum: 2020-03-09 21:13
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRRRR2G3P3A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR2G3P3A0_gravload_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR2G3P3A0_gravload_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3RRRRR2G3P3A0_gravload_para_pf_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRRRR2G3P3A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRRRR2G3P3A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRRRR2G3P3A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR2G3P3A0_gravload_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR2G3P3A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:11:49
% EndTime: 2020-03-09 21:11:50
% DurationCPUTime: 0.81s
% Computational Cost: add. (636->152), mult. (1266->271), div. (72->11), fcn. (666->42), ass. (0->124)
t1104 = cos(qJ(3,1));
t1161 = t1104 ^ 2;
t1101 = cos(qJ(3,2));
t1160 = t1101 ^ 2;
t1098 = cos(qJ(3,3));
t1159 = t1098 ^ 2;
t1158 = m(1) * rSges(1,2);
t1157 = m(3) * rSges(3,3);
t1156 = rSges(2,2) * g(3);
t1155 = rSges(3,2) * g(3);
t1100 = cos(qJ(1,3));
t1154 = pkin(1) * t1100;
t1103 = cos(qJ(1,2));
t1153 = pkin(1) * t1103;
t1106 = cos(qJ(1,1));
t1152 = pkin(1) * t1106;
t1091 = sin(qJ(1,3));
t1151 = t1091 * pkin(1);
t1094 = sin(qJ(1,2));
t1150 = t1094 * pkin(1);
t1097 = sin(qJ(1,1));
t1149 = t1097 * pkin(1);
t1086 = legFrame(3,2);
t1063 = sin(t1086);
t1066 = cos(t1086);
t1040 = t1066 * g(1) - t1063 * g(2);
t1148 = rSges(3,1) * t1040;
t1087 = legFrame(2,2);
t1064 = sin(t1087);
t1067 = cos(t1087);
t1041 = t1067 * g(1) - t1064 * g(2);
t1147 = rSges(3,1) * t1041;
t1088 = legFrame(1,2);
t1065 = sin(t1088);
t1068 = cos(t1088);
t1042 = t1068 * g(1) - t1065 * g(2);
t1146 = rSges(3,1) * t1042;
t1085 = qJ(1,1) + qJ(2,1);
t1084 = qJ(1,2) + qJ(2,2);
t1083 = qJ(1,3) + qJ(2,3);
t1051 = sin(t1083);
t1054 = cos(t1083);
t1075 = 0.1e1 / t1098;
t1089 = sin(qJ(3,3));
t1113 = t1098 * rSges(3,1) - t1089 * rSges(3,2);
t1145 = (-(t1063 * g(1) + t1066 * g(2)) * t1113 + (-g(3) * t1051 + t1040 * t1054) * (t1089 * rSges(3,1) + t1098 * rSges(3,2))) * t1075;
t1052 = sin(t1084);
t1055 = cos(t1084);
t1078 = 0.1e1 / t1101;
t1092 = sin(qJ(3,2));
t1112 = t1101 * rSges(3,1) - t1092 * rSges(3,2);
t1144 = (-(t1064 * g(1) + t1067 * g(2)) * t1112 + (-g(3) * t1052 + t1041 * t1055) * (t1092 * rSges(3,1) + t1101 * rSges(3,2))) * t1078;
t1053 = sin(t1085);
t1056 = cos(t1085);
t1081 = 0.1e1 / t1104;
t1095 = sin(qJ(3,1));
t1111 = t1104 * rSges(3,1) - t1095 * rSges(3,2);
t1143 = (-(t1065 * g(1) + t1068 * g(2)) * t1111 + (-g(3) * t1053 + t1042 * t1056) * (t1095 * rSges(3,1) + t1104 * rSges(3,2))) * t1081;
t1090 = sin(qJ(2,3));
t1099 = cos(qJ(2,3));
t1037 = t1091 * t1090 - t1100 * t1099;
t1142 = t1037 * t1098;
t1093 = sin(qJ(2,2));
t1102 = cos(qJ(2,2));
t1038 = t1094 * t1093 - t1103 * t1102;
t1141 = t1038 * t1101;
t1096 = sin(qJ(2,1));
t1105 = cos(qJ(2,1));
t1039 = t1097 * t1096 - t1106 * t1105;
t1140 = t1039 * t1104;
t1139 = t1063 * t1089;
t1138 = t1064 * t1092;
t1137 = t1065 * t1095;
t1136 = t1066 * t1089;
t1135 = t1067 * t1092;
t1134 = t1068 * t1095;
t1071 = 0.1e1 / t1090;
t1133 = t1071 * t1075;
t1072 = 0.1e1 / t1093;
t1132 = t1072 * t1078;
t1073 = 0.1e1 / t1096;
t1131 = t1073 * t1081;
t1130 = g(3) * t1158;
t1050 = m(1) * rSges(1,1) + m(2) * pkin(1);
t1129 = (m(3) * pkin(1) + t1050) * g(3);
t1128 = pkin(2) * t1037 * t1159;
t1127 = pkin(2) * t1038 * t1160;
t1126 = pkin(2) * t1039 * t1161;
t1125 = t1099 * t1089 * pkin(1);
t1124 = t1102 * t1092 * pkin(1);
t1123 = t1105 * t1095 * pkin(1);
t1108 = rSges(2,1) * g(3);
t1031 = m(2) * (rSges(2,2) * t1040 + t1108);
t1034 = rSges(3,2) * t1040;
t1057 = qJ(3,3) + t1083;
t1058 = -qJ(3,3) + t1083;
t1107 = rSges(3,1) * g(3);
t1069 = g(3) * t1157;
t1116 = (m(2) * (rSges(2,1) * t1040 - t1156) + t1069) * t1051;
t1019 = t1091 * (t1050 * t1040 - t1130) + (-(t1034 - t1107) * cos(t1058) / 0.2e1 - (-t1148 - t1155) * sin(t1058) / 0.2e1 + (t1034 + t1107) * cos(t1057) / 0.2e1 + (t1148 - t1155) * sin(t1057) / 0.2e1 + t1040 * t1151) * m(3) + (-t1040 * t1157 + t1031) * t1054 + t1116 + (t1040 * t1158 + t1129) * t1100;
t1122 = t1019 * t1133;
t1032 = m(2) * (rSges(2,2) * t1041 + t1108);
t1035 = rSges(3,2) * t1041;
t1059 = qJ(3,2) + t1084;
t1060 = -qJ(3,2) + t1084;
t1115 = (m(2) * (rSges(2,1) * t1041 - t1156) + t1069) * t1052;
t1020 = t1094 * (t1050 * t1041 - t1130) + (-(t1035 - t1107) * cos(t1060) / 0.2e1 - (-t1147 - t1155) * sin(t1060) / 0.2e1 + (t1035 + t1107) * cos(t1059) / 0.2e1 + (t1147 - t1155) * sin(t1059) / 0.2e1 + t1041 * t1150) * m(3) + (-t1041 * t1157 + t1032) * t1055 + t1115 + (t1041 * t1158 + t1129) * t1103;
t1121 = t1020 * t1132;
t1033 = m(2) * (rSges(2,2) * t1042 + t1108);
t1036 = rSges(3,2) * t1042;
t1061 = qJ(3,1) + t1085;
t1062 = -qJ(3,1) + t1085;
t1114 = (m(2) * (rSges(2,1) * t1042 - t1156) + t1069) * t1053;
t1021 = t1097 * (t1050 * t1042 - t1130) + (-(t1036 - t1107) * cos(t1062) / 0.2e1 - (-t1146 - t1155) * sin(t1062) / 0.2e1 + (t1036 + t1107) * cos(t1061) / 0.2e1 + (t1146 - t1155) * sin(t1061) / 0.2e1 + t1042 * t1149) * m(3) + (-t1042 * t1157 + t1033) * t1056 + t1114 + (t1042 * t1158 + t1129) * t1106;
t1120 = t1021 * t1131;
t1022 = t1031 * t1054 + t1116 + (t1113 * t1054 * g(3) + (-rSges(3,3) * t1054 + t1051 * t1113) * t1040) * m(3);
t1119 = t1022 * t1071 / t1159;
t1023 = t1032 * t1055 + t1115 + (t1112 * t1055 * g(3) + (-rSges(3,3) * t1055 + t1052 * t1112) * t1041) * m(3);
t1118 = t1023 * t1072 / t1160;
t1024 = t1033 * t1056 + t1114 + (t1111 * t1056 * g(3) + (-rSges(3,3) * t1056 + t1053 * t1111) * t1042) * m(3);
t1117 = t1024 * t1073 / t1161;
t1110 = 0.1e1 / pkin(1);
t1109 = 0.1e1 / pkin(2);
t1 = [-m(4) * g(1) + ((-t1068 * t1140 + t1137) * t1120 + (-t1067 * t1141 + t1138) * t1121 + (-t1066 * t1142 + t1139) * t1122) * t1110 + (((t1068 * t1126 + (-pkin(2) * t1137 - t1068 * t1152) * t1104 - t1065 * t1123) * t1117 + (t1067 * t1127 + (-pkin(2) * t1138 - t1067 * t1153) * t1101 - t1064 * t1124) * t1118 + (t1066 * t1128 + (-pkin(2) * t1139 - t1066 * t1154) * t1098 - t1063 * t1125) * t1119) * t1110 + (t1063 * t1145 + t1064 * t1144 + t1065 * t1143) * m(3)) * t1109; -m(4) * g(2) + ((t1065 * t1140 + t1134) * t1120 + (t1064 * t1141 + t1135) * t1121 + (t1063 * t1142 + t1136) * t1122) * t1110 + (((-t1065 * t1126 + (-pkin(2) * t1134 + t1065 * t1152) * t1104 - t1068 * t1123) * t1117 + (-t1064 * t1127 + (-pkin(2) * t1135 + t1064 * t1153) * t1101 - t1067 * t1124) * t1118 + (-t1063 * t1128 + (-pkin(2) * t1136 + t1063 * t1154) * t1098 - t1066 * t1125) * t1119) * t1110 + (t1066 * t1145 + t1067 * t1144 + t1068 * t1143) * m(3)) * t1109; -m(4) * g(3) + (-t1051 * t1071 * t1019 - t1052 * t1072 * t1020 - t1053 * t1073 * t1021 + ((pkin(2) * (t1106 * t1096 + t1097 * t1105) * t1104 + t1149) * t1024 * t1131 + (pkin(2) * (t1103 * t1093 + t1094 * t1102) * t1101 + t1150) * t1023 * t1132 + (pkin(2) * (t1100 * t1090 + t1091 * t1099) * t1098 + t1151) * t1022 * t1133) * t1109) * t1110;];
taugX  = t1;
