% Calculate minimal parameter regressor of inertia matrix for parallel robot
% P3PRRR2G2A0
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [8x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3PRRR2G2A0_convert_par2_MPV_fixb.m

% Output:
% MMX [3x3]
%   minimal parameter regressor of inertia matrix for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:21
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMX = P3PRRR2G2A0_inertia_para_pf_mdp(xP, qJ, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(2,1),zeros(8,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR2G2A0_inertia_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR2G2A0_inertia_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRR2G2A0_inertia_para_pf_mdp: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR2G2A0_inertia_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR2G2A0_inertia_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [8 1]), ...
  'P3PRRR2G2A0_inertia_para_pf_mdp: MDP has to be [8x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_MMreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:21:44
% EndTime: 2020-03-09 21:21:45
% DurationCPUTime: 1.39s
% Computational Cost: add. (963->209), mult. (1756->413), div. (918->9), fcn. (2226->21), ass. (0->175)
t1228 = MDP(2) / pkin(1) ^ 2;
t1116 = cos(qJ(3,3));
t1117 = cos(qJ(2,3));
t1110 = sin(qJ(3,3));
t1111 = sin(qJ(2,3));
t1188 = t1111 * t1110;
t1071 = (pkin(2) * t1116 + pkin(1)) * t1117 - pkin(2) * t1188;
t1080 = t1117 * t1116 - t1188;
t1227 = t1071 * t1080;
t1122 = 0.1e1 / pkin(2);
t1226 = t1071 * t1122;
t1118 = cos(qJ(3,2));
t1119 = cos(qJ(2,2));
t1112 = sin(qJ(3,2));
t1113 = sin(qJ(2,2));
t1187 = t1113 * t1112;
t1072 = (pkin(2) * t1118 + pkin(1)) * t1119 - pkin(2) * t1187;
t1081 = t1119 * t1118 - t1187;
t1225 = t1072 * t1081;
t1224 = t1072 * t1122;
t1120 = cos(qJ(3,1));
t1121 = cos(qJ(2,1));
t1114 = sin(qJ(3,1));
t1115 = sin(qJ(2,1));
t1186 = t1115 * t1114;
t1073 = (pkin(2) * t1120 + pkin(1)) * t1121 - pkin(2) * t1186;
t1082 = t1121 * t1120 - t1186;
t1223 = t1073 * t1082;
t1222 = t1073 * t1122;
t1102 = 0.1e1 / t1110 ^ 2;
t1221 = t1080 ^ 2 * t1102;
t1104 = 0.1e1 / t1112 ^ 2;
t1220 = t1081 ^ 2 * t1104;
t1106 = 0.1e1 / t1114 ^ 2;
t1219 = t1082 ^ 2 * t1106;
t1107 = legFrame(3,2);
t1095 = sin(t1107);
t1218 = t1080 * t1095;
t1098 = cos(t1107);
t1217 = t1080 * t1098;
t1108 = legFrame(2,2);
t1096 = sin(t1108);
t1216 = t1081 * t1096;
t1099 = cos(t1108);
t1215 = t1081 * t1099;
t1109 = legFrame(1,2);
t1097 = sin(t1109);
t1214 = t1082 * t1097;
t1100 = cos(t1109);
t1213 = t1082 * t1100;
t1092 = sin(qJ(2,3) + qJ(3,3));
t1083 = t1111 * pkin(1) + pkin(2) * t1092;
t1101 = 0.1e1 / t1110;
t1212 = t1083 * t1101;
t1093 = sin(qJ(2,2) + qJ(3,2));
t1084 = t1113 * pkin(1) + pkin(2) * t1093;
t1103 = 0.1e1 / t1112;
t1211 = t1084 * t1103;
t1094 = sin(qJ(2,1) + qJ(3,1));
t1085 = t1115 * pkin(1) + pkin(2) * t1094;
t1105 = 0.1e1 / t1114;
t1210 = t1085 * t1105;
t1209 = t1092 * t1101;
t1208 = t1093 * t1103;
t1207 = t1094 * t1105;
t1206 = t1095 * t1101;
t1205 = t1096 * t1103;
t1204 = t1097 * t1105;
t1203 = t1098 * t1095;
t1202 = t1098 * t1101;
t1201 = t1099 * t1096;
t1200 = t1099 * t1103;
t1199 = t1100 * t1097;
t1198 = t1100 * t1105;
t1197 = t1101 * t1116;
t1123 = 0.1e1 / pkin(1);
t1196 = t1101 * t1123;
t1195 = t1102 * t1116;
t1194 = t1103 * t1118;
t1193 = t1103 * t1123;
t1192 = t1104 * t1118;
t1191 = t1105 * t1120;
t1190 = t1105 * t1123;
t1189 = t1106 * t1120;
t1185 = t1122 * t1123;
t1184 = (t1199 * t1219 + t1201 * t1220 + t1203 * t1221) * t1228 + (-t1199 - t1201 - t1203) * MDP(1);
t1183 = t1101 * t1227;
t1182 = t1071 * t1206;
t1181 = t1071 * t1202;
t1180 = t1103 * t1225;
t1179 = t1072 * t1205;
t1178 = t1072 * t1200;
t1177 = t1105 * t1223;
t1176 = t1073 * t1204;
t1175 = t1073 * t1198;
t1174 = t1080 * t1092 * t1102;
t1173 = t1080 * t1206;
t1172 = t1080 * t1202;
t1171 = t1080 * t1197;
t1170 = t1080 * t1195;
t1169 = t1081 * t1093 * t1104;
t1168 = t1081 * t1205;
t1167 = t1081 * t1200;
t1166 = t1081 * t1194;
t1165 = t1081 * t1192;
t1164 = t1082 * t1094 * t1106;
t1163 = t1082 * t1204;
t1162 = t1082 * t1198;
t1161 = t1082 * t1191;
t1160 = t1082 * t1189;
t1159 = t1092 * t1197;
t1158 = t1092 * t1195;
t1157 = t1093 * t1194;
t1156 = t1093 * t1192;
t1155 = t1094 * t1191;
t1154 = t1094 * t1189;
t1153 = t1101 * t1185;
t1152 = t1103 * t1185;
t1151 = t1105 * t1185;
t1150 = t1092 * t1196;
t1149 = t1093 * t1193;
t1148 = t1094 * t1190;
t1147 = t1151 * t1199 * t1223 + t1152 * t1201 * t1225 + t1153 * t1203 * t1227;
t1146 = t1071 * t1170;
t1145 = t1071 * t1158;
t1144 = t1072 * t1165;
t1143 = t1072 * t1156;
t1142 = t1073 * t1160;
t1141 = t1073 * t1154;
t1140 = t1095 * t1171;
t1139 = t1098 * t1171;
t1138 = t1098 * t1170;
t1137 = t1096 * t1166;
t1136 = t1099 * t1166;
t1135 = t1099 * t1165;
t1134 = t1097 * t1161;
t1133 = t1100 * t1161;
t1132 = t1100 * t1160;
t1131 = (t1080 - t1226) * t1196;
t1130 = (0.2e1 * t1080 - t1226) * t1196;
t1129 = (t1081 - t1224) * t1193;
t1128 = (0.2e1 * t1081 - t1224) * t1193;
t1127 = (t1082 - t1222) * t1190;
t1126 = (0.2e1 * t1082 - t1222) * t1190;
t1125 = (-t1071 * t1095 * t1138 - t1072 * t1096 * t1135 - t1073 * t1097 * t1132) * MDP(6);
t1091 = t1100 ^ 2;
t1090 = t1099 ^ 2;
t1089 = t1098 ^ 2;
t1088 = t1097 ^ 2;
t1087 = t1096 ^ 2;
t1086 = t1095 ^ 2;
t1076 = t1085 * t1151;
t1075 = t1084 * t1152;
t1074 = t1083 * t1153;
t1070 = t1076 - t1148;
t1069 = t1075 - t1149;
t1068 = t1074 - t1150;
t1067 = t1076 - 0.2e1 * t1148;
t1066 = t1075 - 0.2e1 * t1149;
t1065 = t1074 - 0.2e1 * t1150;
t1060 = t1100 * t1127;
t1059 = t1099 * t1129;
t1058 = t1098 * t1131;
t1057 = t1097 * t1127;
t1056 = t1096 * t1129;
t1055 = t1095 * t1131;
t1054 = t1100 * t1126;
t1053 = t1099 * t1128;
t1052 = t1098 * t1130;
t1051 = t1097 * t1126;
t1050 = t1096 * t1128;
t1049 = t1095 * t1130;
t1047 = (-t1098 * t1174 - t1099 * t1169 - t1100 * t1164) * t1228;
t1046 = (-t1095 * t1174 - t1096 * t1169 - t1097 * t1164) * t1228;
t1 = [(t1091 + t1090 + t1089) * MDP(1) + (t1049 * t1140 + t1050 * t1137 + t1051 * t1134) * MDP(6) + (-t1049 * t1218 - t1050 * t1216 - t1051 * t1214) * MDP(7) + MDP(8) + (t1086 * t1221 + t1087 * t1220 + t1088 * t1219) * t1228 + ((t1055 * t1173 + t1056 * t1168 + t1057 * t1163) * MDP(5) + ((-t1055 * t1182 - t1056 * t1179 - t1057 * t1176) * MDP(5) + (-t1086 * t1146 - t1087 * t1144 - t1088 * t1142) * MDP(6) + (t1086 * t1183 + t1087 * t1180 + t1088 * t1177) * MDP(7)) * t1122) * t1123; (t1052 * t1140 + t1053 * t1137 + t1054 * t1134) * MDP(6) + (-t1052 * t1218 - t1053 * t1216 - t1054 * t1214 + t1147) * MDP(7) + ((t1058 * t1173 + t1059 * t1168 + t1060 * t1163) * MDP(5) + ((-t1058 * t1182 - t1059 * t1179 - t1060 * t1176) * MDP(5) + t1125) * t1122) * t1123 + t1184; t1046 + (t1065 * t1140 + t1066 * t1137 + t1067 * t1134) * MDP(6) + (-t1065 * t1218 - t1066 * t1216 - t1067 * t1214) * MDP(7) + ((t1068 * t1173 + t1069 * t1168 + t1070 * t1163) * MDP(5) + ((-t1068 * t1182 - t1069 * t1179 - t1070 * t1176) * MDP(5) + (t1095 * t1145 + t1096 * t1143 + t1097 * t1141) * MDP(6) + (-t1092 * t1182 - t1093 * t1179 - t1094 * t1176) * MDP(7)) * t1122) * t1123; (t1049 * t1139 + t1050 * t1136 + t1051 * t1133) * MDP(6) + (-t1049 * t1217 - t1050 * t1215 - t1051 * t1213 + t1147) * MDP(7) + ((t1055 * t1172 + t1056 * t1167 + t1057 * t1162) * MDP(5) + ((-t1055 * t1181 - t1056 * t1178 - t1057 * t1175) * MDP(5) + t1125) * t1122) * t1123 + t1184; (t1088 + t1087 + t1086) * MDP(1) + (t1052 * t1139 + t1053 * t1136 + t1054 * t1133) * MDP(6) + (-t1052 * t1217 - t1053 * t1215 - t1054 * t1213) * MDP(7) + MDP(8) + (t1089 * t1221 + t1090 * t1220 + t1091 * t1219) * t1228 + ((t1058 * t1172 + t1059 * t1167 + t1060 * t1162) * MDP(5) + ((-t1058 * t1181 - t1059 * t1178 - t1060 * t1175) * MDP(5) + (-t1089 * t1146 - t1090 * t1144 - t1091 * t1142) * MDP(6) + (t1089 * t1183 + t1090 * t1180 + t1091 * t1177) * MDP(7)) * t1122) * t1123; t1047 + (t1065 * t1139 + t1066 * t1136 + t1067 * t1133) * MDP(6) + (-t1065 * t1217 - t1066 * t1215 - t1067 * t1213) * MDP(7) + ((t1068 * t1172 + t1069 * t1167 + t1070 * t1162) * MDP(5) + ((-t1068 * t1181 - t1069 * t1178 - t1070 * t1175) * MDP(5) + (t1098 * t1145 + t1099 * t1143 + t1100 * t1141) * MDP(6) + (-t1092 * t1181 - t1093 * t1178 - t1094 * t1175) * MDP(7)) * t1122) * t1123; t1046 + (-t1049 * t1159 - t1050 * t1157 - t1051 * t1155) * MDP(6) + (t1092 * t1049 + t1093 * t1050 + t1094 * t1051) * MDP(7) + ((-t1055 * t1209 - t1056 * t1208 - t1057 * t1207) * MDP(5) + ((t1055 * t1212 + t1056 * t1211 + t1057 * t1210) * MDP(5) + (t1083 * t1095 * t1170 + t1084 * t1096 * t1165 + t1085 * t1097 * t1160) * MDP(6) + (-t1083 * t1173 - t1084 * t1168 - t1085 * t1163) * MDP(7)) * t1122) * t1123; t1047 + (-t1052 * t1159 - t1053 * t1157 - t1054 * t1155) * MDP(6) + (t1092 * t1052 + t1093 * t1053 + t1094 * t1054) * MDP(7) + ((-t1058 * t1209 - t1059 * t1208 - t1060 * t1207) * MDP(5) + ((t1058 * t1212 + t1059 * t1211 + t1060 * t1210) * MDP(5) + (t1083 * t1138 + t1084 * t1135 + t1085 * t1132) * MDP(6) + (-t1083 * t1172 - t1084 * t1167 - t1085 * t1162) * MDP(7)) * t1122) * t1123; (-t1065 * t1159 - t1066 * t1157 - t1067 * t1155) * MDP(6) + (t1092 * t1065 + t1093 * t1066 + t1094 * t1067) * MDP(7) + MDP(8) + (t1092 ^ 2 * t1102 + t1093 ^ 2 * t1104 + t1094 ^ 2 * t1106) * t1228 + ((-t1068 * t1209 - t1069 * t1208 - t1070 * t1207) * MDP(5) + ((t1068 * t1212 + t1069 * t1211 + t1070 * t1210) * MDP(5) + (-t1083 * t1158 - t1084 * t1156 - t1085 * t1154) * MDP(6) + (t1083 * t1209 + t1084 * t1208 + t1085 * t1207) * MDP(7)) * t1122) * t1123;];
%% Postprocessing: Reshape Output
% From vec2mat_3_matlab.m
res = [t1(1), t1(2), t1(3); t1(4), t1(5), t1(6); t1(7), t1(8), t1(9);];
MMX  = res;
