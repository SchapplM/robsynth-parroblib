% Calculate minimal parameter regressor of vector of centrifugal and coriolis load for parallel robot
% P3RPR1G1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% qJ [2x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d3]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [10x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RPR1G1A0_convert_par2_MPV_fixb.m

% Output:
% taucX [3x1]
%   minimal parameter regressor of vector of coriolis and centrifugal joint torques
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 14:58
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RPR1G1A0_coriolisvec_para_pf_mdp(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(2,3),zeros(3,3),zeros(3,3),zeros(4,1),zeros(10,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPR1G1A0_coriolisvec_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPR1G1A0_coriolisvec_para_pf_mdp: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3RPR1G1A0_coriolisvec_para_pf_mdp: qJ has to be [2x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RPR1G1A0_coriolisvec_para_pf_mdp: pkin has to be [4x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPR1G1A0_coriolisvec_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPR1G1A0_coriolisvec_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'P3RPR1G1A0_coriolisvec_para_pf_mdp: MDP has to be [10x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_tauCreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:58:17
% EndTime: 2019-05-03 14:58:18
% DurationCPUTime: 1.80s
% Computational Cost: add. (11107->205), mult. (17565->368), div. (567->9), fcn. (8218->14), ass. (0->153)
t1150 = pkin(1) + pkin(2);
t1164 = koppelP(3,2);
t1167 = koppelP(3,1);
t1120 = -qJ(2,3) * t1164 + t1150 * t1167;
t1123 = qJ(2,3) * t1167 + t1150 * t1164;
t1148 = xDP(2);
t1126 = t1148 * t1150;
t1151 = xP(3);
t1134 = sin(t1151);
t1135 = cos(t1151);
t1147 = xDP(3);
t1149 = xDP(1);
t1072 = qJ(2,3) * t1149 + t1126 + (t1120 * t1135 - t1123 * t1134) * t1147;
t1127 = t1149 * t1150;
t1075 = qJ(2,3) * t1148 - t1127 + (t1120 * t1134 + t1123 * t1135) * t1147;
t1138 = legFrame(3,3);
t1128 = sin(t1138);
t1131 = cos(t1138);
t1141 = sin(qJ(1,3));
t1144 = cos(qJ(1,3));
t1057 = (t1072 * t1141 - t1075 * t1144) * t1131 + (t1072 * t1144 + t1075 * t1141) * t1128;
t1102 = t1134 * t1167 + t1135 * t1164;
t1090 = t1102 * t1147 - t1149;
t1105 = -t1134 * t1164 + t1135 * t1167;
t1093 = t1105 * t1147 + t1148;
t1066 = (-t1090 * t1144 + t1093 * t1141) * t1131 + t1128 * (t1090 * t1141 + t1093 * t1144);
t1114 = -qJ(2,3) * t1144 + t1141 * t1150;
t1117 = qJ(2,3) * t1141 + t1144 * t1150;
t1078 = t1114 * t1131 + t1117 * t1128;
t1081 = -t1114 * t1128 + t1117 * t1131;
t1136 = t1147 ^ 2;
t1152 = qJ(2,3) ^ 2;
t1153 = 0.1e1 / qJ(2,3);
t1154 = 0.1e1 / qJ(2,3) ^ 2;
t1171 = (pkin(1) ^ 2);
t1187 = -t1171 + (-2 * pkin(1) - pkin(2)) * pkin(2);
t1155 = t1153 / t1152;
t1215 = t1057 * t1155;
t1180 = -t1066 * (((-t1152 + t1187) * t1066 + t1057 * t1150) * t1153 * t1154 + t1150 * t1215) - (t1078 * t1102 + t1081 * t1105) * t1153 * t1136;
t1165 = koppelP(2,2);
t1168 = koppelP(2,1);
t1121 = -qJ(2,2) * t1165 + t1150 * t1168;
t1124 = qJ(2,2) * t1168 + t1150 * t1165;
t1073 = qJ(2,2) * t1149 + t1126 + (t1121 * t1135 - t1124 * t1134) * t1147;
t1076 = qJ(2,2) * t1148 - t1127 + (t1121 * t1134 + t1124 * t1135) * t1147;
t1139 = legFrame(2,3);
t1129 = sin(t1139);
t1132 = cos(t1139);
t1142 = sin(qJ(1,2));
t1145 = cos(qJ(1,2));
t1058 = (t1073 * t1142 - t1076 * t1145) * t1132 + (t1073 * t1145 + t1076 * t1142) * t1129;
t1103 = t1134 * t1168 + t1135 * t1165;
t1091 = t1103 * t1147 - t1149;
t1106 = -t1134 * t1165 + t1135 * t1168;
t1094 = t1106 * t1147 + t1148;
t1067 = (-t1091 * t1145 + t1094 * t1142) * t1132 + t1129 * (t1091 * t1142 + t1094 * t1145);
t1115 = -qJ(2,2) * t1145 + t1142 * t1150;
t1118 = qJ(2,2) * t1142 + t1145 * t1150;
t1079 = t1115 * t1132 + t1118 * t1129;
t1082 = -t1115 * t1129 + t1118 * t1132;
t1156 = qJ(2,2) ^ 2;
t1157 = 0.1e1 / qJ(2,2);
t1158 = 0.1e1 / qJ(2,2) ^ 2;
t1159 = t1157 / t1156;
t1214 = t1058 * t1159;
t1179 = -t1067 * (((-t1156 + t1187) * t1067 + t1058 * t1150) * t1157 * t1158 + t1150 * t1214) - (t1079 * t1103 + t1082 * t1106) * t1157 * t1136;
t1166 = koppelP(1,2);
t1169 = koppelP(1,1);
t1122 = -qJ(2,1) * t1166 + t1150 * t1169;
t1125 = qJ(2,1) * t1169 + t1150 * t1166;
t1074 = qJ(2,1) * t1149 + t1126 + (t1122 * t1135 - t1125 * t1134) * t1147;
t1077 = qJ(2,1) * t1148 - t1127 + (t1122 * t1134 + t1125 * t1135) * t1147;
t1140 = legFrame(1,3);
t1130 = sin(t1140);
t1133 = cos(t1140);
t1143 = sin(qJ(1,1));
t1146 = cos(qJ(1,1));
t1059 = (t1074 * t1143 - t1077 * t1146) * t1133 + (t1074 * t1146 + t1077 * t1143) * t1130;
t1104 = t1134 * t1169 + t1135 * t1166;
t1092 = t1104 * t1147 - t1149;
t1107 = -t1134 * t1166 + t1135 * t1169;
t1095 = t1107 * t1147 + t1148;
t1068 = (-t1092 * t1146 + t1095 * t1143) * t1133 + t1130 * (t1092 * t1143 + t1095 * t1146);
t1116 = -qJ(2,1) * t1146 + t1143 * t1150;
t1119 = qJ(2,1) * t1143 + t1146 * t1150;
t1080 = t1116 * t1133 + t1119 * t1130;
t1083 = -t1116 * t1130 + t1119 * t1133;
t1160 = qJ(2,1) ^ 2;
t1161 = 0.1e1 / qJ(2,1);
t1162 = 0.1e1 / qJ(2,1) ^ 2;
t1163 = t1161 / t1160;
t1213 = t1059 * t1163;
t1178 = -t1068 * (((-t1160 + t1187) * t1068 + t1059 * t1150) * t1161 * t1162 + t1150 * t1213) - (t1080 * t1104 + t1083 * t1107) * t1161 * t1136;
t1219 = 2 * MDP(6);
t1096 = t1128 * t1144 + t1131 * t1141;
t1097 = -t1128 * t1141 + t1131 * t1144;
t1190 = t1066 * t1215;
t1203 = t1066 * t1154;
t1051 = -t1190 + (-(-t1066 * t1150 + t1057) * t1203 + (-t1096 * t1102 - t1097 * t1105) * t1136) * t1153;
t1218 = pkin(1) * t1051;
t1098 = t1129 * t1145 + t1132 * t1142;
t1099 = -t1129 * t1142 + t1132 * t1145;
t1189 = t1067 * t1214;
t1202 = t1067 * t1158;
t1052 = -t1189 + (-(-t1067 * t1150 + t1058) * t1202 + (-t1098 * t1103 - t1099 * t1106) * t1136) * t1157;
t1217 = pkin(1) * t1052;
t1100 = t1130 * t1146 + t1133 * t1143;
t1101 = -t1130 * t1143 + t1133 * t1146;
t1188 = t1068 * t1213;
t1201 = t1068 * t1162;
t1053 = -t1188 + (-(-t1068 * t1150 + t1059) * t1201 + (-t1100 * t1104 - t1101 * t1107) * t1136) * t1161;
t1216 = pkin(1) * t1053;
t1084 = t1120 * t1144 + t1123 * t1141;
t1087 = t1120 * t1141 - t1123 * t1144;
t1060 = (-t1084 * t1134 + t1087 * t1135) * t1131 + (t1084 * t1135 + t1087 * t1134) * t1128;
t1212 = t1060 * t1153;
t1085 = t1121 * t1145 + t1124 * t1142;
t1088 = t1121 * t1142 - t1124 * t1145;
t1061 = (-t1085 * t1134 + t1088 * t1135) * t1132 + (t1085 * t1135 + t1088 * t1134) * t1129;
t1211 = t1061 * t1157;
t1086 = t1122 * t1146 + t1125 * t1143;
t1089 = t1122 * t1143 - t1125 * t1146;
t1062 = (-t1086 * t1134 + t1089 * t1135) * t1133 + (t1086 * t1135 + t1089 * t1134) * t1130;
t1210 = t1062 * t1161;
t1063 = t1066 ^ 2;
t1209 = t1063 * t1154;
t1208 = t1063 * t1155;
t1064 = t1067 ^ 2;
t1207 = t1064 * t1158;
t1206 = t1064 * t1159;
t1065 = t1068 ^ 2;
t1205 = t1065 * t1162;
t1204 = t1065 * t1163;
t1200 = t1078 * t1153;
t1199 = t1079 * t1157;
t1198 = t1080 * t1161;
t1197 = t1081 * t1153;
t1196 = t1082 * t1157;
t1195 = t1083 * t1161;
t1191 = 2 * MDP(5);
t1177 = (t1051 + t1190) * t1191 + t1057 * t1203 * t1219 + (MDP(1) * t1051 + MDP(4) * (-t1180 + 0.2e1 * t1218) + MDP(6) * ((t1152 + t1171) * t1051 - pkin(1) * t1180)) * t1153;
t1176 = (t1052 + t1189) * t1191 + t1058 * t1202 * t1219 + (MDP(1) * t1052 + MDP(4) * (-t1179 + 0.2e1 * t1217) + MDP(6) * ((t1156 + t1171) * t1052 - pkin(1) * t1179)) * t1157;
t1175 = (t1053 + t1188) * t1191 + t1059 * t1201 * t1219 + (MDP(1) * t1053 + MDP(4) * (-t1178 + 0.2e1 * t1216) + MDP(6) * ((t1160 + t1171) * t1053 - pkin(1) * t1178)) * t1161;
t1113 = t1143 * t1166 + t1146 * t1169;
t1112 = t1143 * t1169 - t1146 * t1166;
t1111 = t1142 * t1165 + t1145 * t1168;
t1110 = t1142 * t1168 - t1145 * t1165;
t1109 = t1141 * t1164 + t1144 * t1167;
t1108 = t1141 * t1167 - t1144 * t1164;
t1050 = t1178 - t1216;
t1048 = t1179 - t1217;
t1046 = t1180 - t1218;
t1 = [(-t1051 * t1197 - t1052 * t1196 - t1053 * t1195) * MDP(4) + (-t1081 * t1208 - t1082 * t1206 - t1083 * t1204) * MDP(5) + (t1046 * t1197 + t1048 * t1196 + t1050 * t1195 - t1081 * t1209 - t1082 * t1207 - t1083 * t1205) * MDP(6) + (-MDP(8) * t1135 + MDP(9) * t1134) * t1136 + t1175 * t1101 + t1176 * t1099 + t1177 * t1097; (-t1051 * t1200 - t1052 * t1199 - t1053 * t1198) * MDP(4) + (-t1078 * t1208 - t1079 * t1206 - t1080 * t1204) * MDP(5) + (t1046 * t1200 + t1048 * t1199 + t1050 * t1198 - t1078 * t1209 - t1079 * t1207 - t1080 * t1205) * MDP(6) + (-MDP(8) * t1134 - MDP(9) * t1135) * t1136 + t1175 * t1100 + t1176 * t1098 + t1177 * t1096; (-t1051 * t1212 - t1052 * t1211 - t1053 * t1210) * MDP(4) + (-t1060 * t1208 - t1061 * t1206 - t1062 * t1204) * MDP(5) + (t1046 * t1212 + t1048 * t1211 + t1050 * t1210 - t1060 * t1209 - t1061 * t1207 - t1062 * t1205) * MDP(6) + t1175 * ((t1112 * t1135 - t1113 * t1134) * t1133 + t1130 * (t1112 * t1134 + t1113 * t1135)) + t1176 * ((t1110 * t1135 - t1111 * t1134) * t1132 + t1129 * (t1110 * t1134 + t1111 * t1135)) + t1177 * ((t1108 * t1135 - t1109 * t1134) * t1131 + t1128 * (t1108 * t1134 + t1109 * t1135));];
taucX  = t1;
