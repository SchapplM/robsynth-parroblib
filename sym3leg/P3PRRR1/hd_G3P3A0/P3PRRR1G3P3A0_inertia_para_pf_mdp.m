% Calculate minimal parameter regressor of inertia matrix for parallel robot
% P3PRRR1G3P3A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [8x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3PRRR1G3P3A0_convert_par2_MPV_fixb.m

% Output:
% MMX [3x3]
%   minimal parameter regressor of inertia matrix for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:07
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMX = P3PRRR1G3P3A0_inertia_para_pf_mdp(xP, qJ, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(8,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR1G3P3A0_inertia_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR1G3P3A0_inertia_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3PRRR1G3P3A0_inertia_para_pf_mdp: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR1G3P3A0_inertia_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR1G3P3A0_inertia_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [8 1]), ...
  'P3PRRR1G3P3A0_inertia_para_pf_mdp: MDP has to be [8x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_MMreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:07:01
% EndTime: 2020-03-09 21:07:03
% DurationCPUTime: 2.21s
% Computational Cost: add. (4710->224), mult. (2059->448), div. (981->9), fcn. (2856->24), ass. (0->211)
t1178 = legFrame(3,2);
t1169 = sin(t1178);
t1175 = pkin(7) + qJ(2,3);
t1166 = qJ(3,3) + t1175;
t1148 = sin(t1166);
t1151 = cos(t1166);
t1154 = sin(t1175);
t1157 = cos(t1175);
t1136 = t1148 * t1157 - t1151 * t1154;
t1314 = 0.1e1 / t1136;
t1320 = t1169 * t1314;
t1179 = legFrame(2,2);
t1170 = sin(t1179);
t1176 = pkin(7) + qJ(2,2);
t1167 = qJ(3,2) + t1176;
t1149 = sin(t1167);
t1152 = cos(t1167);
t1155 = sin(t1176);
t1158 = cos(t1176);
t1137 = t1149 * t1158 - t1152 * t1155;
t1313 = 0.1e1 / t1137;
t1319 = t1170 * t1313;
t1180 = legFrame(1,2);
t1171 = sin(t1180);
t1177 = pkin(7) + qJ(2,1);
t1168 = qJ(3,1) + t1177;
t1150 = sin(t1168);
t1153 = cos(t1168);
t1156 = sin(t1177);
t1159 = cos(t1177);
t1138 = t1150 * t1159 - t1153 * t1156;
t1312 = 0.1e1 / t1138;
t1318 = t1171 * t1312;
t1317 = t1312 * t1150;
t1316 = t1313 * t1149;
t1315 = t1314 * t1148;
t1311 = MDP(2) / pkin(2) ^ 2;
t1310 = t1314 ^ 2;
t1139 = pkin(2) * t1154 + pkin(3) * t1148;
t1309 = t1314 * t1139;
t1188 = 0.1e1 / pkin(2);
t1307 = t1314 * t1188;
t1128 = 0.1e1 / t1136 ^ 2;
t1306 = t1128 * t1139;
t1305 = t1128 * t1148;
t1304 = t1313 ^ 2;
t1140 = pkin(2) * t1155 + pkin(3) * t1149;
t1303 = t1313 * t1140;
t1301 = t1313 * t1188;
t1130 = 0.1e1 / t1137 ^ 2;
t1300 = t1130 * t1140;
t1299 = t1130 * t1149;
t1298 = t1312 ^ 2;
t1141 = pkin(2) * t1156 + pkin(3) * t1150;
t1297 = t1312 * t1141;
t1295 = t1312 * t1188;
t1132 = 0.1e1 / t1138 ^ 2;
t1294 = t1132 * t1141;
t1293 = t1132 * t1150;
t1172 = cos(t1178);
t1291 = t1314 * t1172;
t1173 = cos(t1179);
t1289 = t1313 * t1173;
t1174 = cos(t1180);
t1287 = t1312 * t1174;
t1142 = pkin(2) * t1157 + pkin(3) * t1151;
t1286 = t1142 * t1151;
t1285 = t1142 * t1169;
t1143 = pkin(2) * t1158 + pkin(3) * t1152;
t1284 = t1143 * t1152;
t1283 = t1143 * t1170;
t1144 = pkin(2) * t1159 + pkin(3) * t1153;
t1282 = t1144 * t1153;
t1281 = t1144 * t1171;
t1280 = t1151 * t1172;
t1279 = t1152 * t1173;
t1278 = t1153 * t1174;
t1277 = t1169 * t1172;
t1276 = t1170 * t1173;
t1275 = t1171 * t1174;
t1187 = 0.1e1 / pkin(3);
t1274 = t1187 * t1188;
t1273 = t1314 * t1285;
t1272 = t1148 * t1307;
t1271 = t1314 * t1280;
t1270 = t1151 * t1307;
t1181 = sin(qJ(3,3));
t1269 = t1181 * t1306;
t1184 = cos(qJ(3,3));
t1268 = t1184 * t1306;
t1267 = t1313 * t1283;
t1266 = t1149 * t1301;
t1265 = t1313 * t1279;
t1264 = t1152 * t1301;
t1182 = sin(qJ(3,2));
t1263 = t1182 * t1300;
t1185 = cos(qJ(3,2));
t1262 = t1185 * t1300;
t1261 = t1312 * t1281;
t1260 = t1150 * t1295;
t1259 = t1312 * t1278;
t1258 = t1153 * t1295;
t1183 = sin(qJ(3,1));
t1257 = t1183 * t1294;
t1186 = cos(qJ(3,1));
t1256 = t1186 * t1294;
t1255 = t1142 * t1291;
t1254 = t1181 * t1315;
t1253 = t1184 * t1315;
t1252 = t1151 * t1320;
t1251 = t1143 * t1289;
t1250 = t1182 * t1316;
t1249 = t1185 * t1316;
t1248 = t1152 * t1319;
t1247 = t1144 * t1287;
t1246 = t1183 * t1317;
t1245 = t1186 * t1317;
t1244 = t1153 * t1318;
t1243 = t1142 * t1274;
t1242 = t1143 * t1274;
t1241 = t1144 * t1274;
t1145 = t1151 ^ 2;
t1240 = t1145 * t1277;
t1146 = t1152 ^ 2;
t1239 = t1146 * t1276;
t1147 = t1153 ^ 2;
t1238 = t1147 * t1275;
t1237 = t1291 * t1315;
t1236 = t1169 * t1270;
t1235 = t1181 * t1271;
t1234 = t1184 * t1271;
t1233 = t1172 * t1270;
t1232 = t1151 * t1269;
t1231 = t1151 * t1268;
t1230 = t1285 * t1305;
t1160 = t1169 ^ 2;
t1229 = t1128 * t1160 * t1286;
t1228 = t1289 * t1316;
t1227 = t1170 * t1264;
t1226 = t1182 * t1265;
t1225 = t1185 * t1265;
t1224 = t1173 * t1264;
t1223 = t1152 * t1263;
t1222 = t1152 * t1262;
t1221 = t1283 * t1299;
t1161 = t1170 ^ 2;
t1220 = t1130 * t1161 * t1284;
t1219 = t1287 * t1317;
t1218 = t1171 * t1258;
t1217 = t1183 * t1259;
t1216 = t1186 * t1259;
t1215 = t1174 * t1258;
t1214 = t1153 * t1257;
t1213 = t1153 * t1256;
t1212 = t1281 * t1293;
t1162 = t1171 ^ 2;
t1211 = t1132 * t1162 * t1282;
t1210 = t1181 * t1252;
t1209 = t1184 * t1252;
t1208 = t1182 * t1248;
t1207 = t1185 * t1248;
t1206 = t1183 * t1244;
t1205 = t1186 * t1244;
t1204 = t1277 * t1286;
t1203 = t1276 * t1284;
t1202 = t1275 * t1282;
t1201 = t1142 * t1237;
t1163 = t1172 ^ 2;
t1200 = t1163 * t1286 * t1310;
t1199 = t1128 * t1204;
t1198 = t1143 * t1228;
t1164 = t1173 ^ 2;
t1197 = t1164 * t1284 * t1304;
t1196 = t1130 * t1203;
t1195 = t1144 * t1219;
t1165 = t1174 ^ 2;
t1194 = t1165 * t1282 * t1298;
t1193 = t1132 * t1202;
t1192 = t1204 * t1310;
t1191 = t1203 * t1304;
t1190 = t1202 * t1298;
t1126 = (t1275 + t1276 + t1277) * MDP(1);
t1125 = t1274 * t1297;
t1124 = t1274 * t1303;
t1123 = t1274 * t1309;
t1122 = t1241 * t1318;
t1121 = t1242 * t1319;
t1120 = t1243 * t1320;
t1119 = t1241 * t1287;
t1118 = t1242 * t1289;
t1117 = t1243 * t1291;
t1116 = t1125 - 0.2e1 * t1260;
t1115 = t1125 - t1260;
t1114 = t1124 - 0.2e1 * t1266;
t1113 = t1124 - t1266;
t1112 = t1123 - 0.2e1 * t1272;
t1111 = t1123 - t1272;
t1110 = -t1119 + t1215;
t1109 = -t1119 + 0.2e1 * t1215;
t1108 = -t1118 + t1224;
t1107 = -t1118 + 0.2e1 * t1224;
t1106 = -t1117 + t1233;
t1105 = -t1117 + 0.2e1 * t1233;
t1104 = t1122 - 0.2e1 * t1218;
t1103 = t1122 - t1218;
t1102 = t1121 - 0.2e1 * t1227;
t1101 = t1121 - t1227;
t1100 = t1120 - 0.2e1 * t1236;
t1099 = t1120 - t1236;
t1098 = (t1244 * t1317 + t1248 * t1316 + t1252 * t1315) * t1311;
t1 = [(t1162 + t1161 + t1160) * MDP(1) + (t1105 * t1234 + t1107 * t1225 + t1109 * t1216) * MDP(6) + (-t1105 * t1235 - t1107 * t1226 - t1109 * t1217) * MDP(7) + MDP(8) + (t1128 * t1145 * t1163 + t1130 * t1146 * t1164 + t1132 * t1147 * t1165) * t1311 + ((t1106 * t1271 + t1108 * t1265 + t1110 * t1259) * MDP(5) + ((-t1106 * t1255 - t1108 * t1251 - t1110 * t1247) * MDP(5) + (-t1184 * t1200 - t1185 * t1197 - t1186 * t1194) * MDP(6) + (t1181 * t1200 + t1182 * t1197 + t1183 * t1194) * MDP(7)) * t1187) * t1188; t1126 + (t1100 * t1234 + t1102 * t1225 + t1104 * t1216) * MDP(6) + (-t1100 * t1235 - t1102 * t1226 - t1104 * t1217) * MDP(7) + (-t1128 * t1240 - t1130 * t1239 - t1132 * t1238) * t1311 + ((t1099 * t1271 + t1101 * t1265 + t1103 * t1259) * MDP(5) + ((-t1099 * t1255 - t1101 * t1251 - t1103 * t1247) * MDP(5) + (t1184 * t1192 + t1185 * t1191 + t1186 * t1190) * MDP(6) + (-t1181 * t1192 - t1182 * t1191 - t1183 * t1190) * MDP(7)) * t1187) * t1188; (t1112 * t1234 + t1114 * t1225 + t1116 * t1216) * MDP(6) + (-t1112 * t1235 - t1114 * t1226 - t1116 * t1217) * MDP(7) + (-t1278 * t1293 - t1279 * t1299 - t1280 * t1305) * t1311 + ((t1111 * t1271 + t1113 * t1265 + t1115 * t1259) * MDP(5) + ((-t1111 * t1255 - t1113 * t1251 - t1115 * t1247) * MDP(5) + (t1184 * t1201 + t1185 * t1198 + t1186 * t1195) * MDP(6) + (-t1181 * t1201 - t1182 * t1198 - t1183 * t1195) * MDP(7)) * t1187) * t1188; t1126 + (-t1105 * t1209 - t1107 * t1207 - t1109 * t1205) * MDP(6) + (t1105 * t1210 + t1107 * t1208 + t1109 * t1206) * MDP(7) + (-t1238 * t1298 - t1239 * t1304 - t1240 * t1310) * t1311 + ((-t1106 * t1252 - t1108 * t1248 - t1110 * t1244) * MDP(5) + ((t1106 * t1273 + t1108 * t1267 + t1110 * t1261) * MDP(5) + (t1184 * t1199 + t1185 * t1196 + t1186 * t1193) * MDP(6) + (-t1181 * t1199 - t1182 * t1196 - t1183 * t1193) * MDP(7)) * t1187) * t1188; (t1165 + t1164 + t1163) * MDP(1) + (-t1100 * t1209 - t1102 * t1207 - t1104 * t1205) * MDP(6) + (t1100 * t1210 + t1102 * t1208 + t1104 * t1206) * MDP(7) + MDP(8) + (t1145 * t1160 * t1310 + t1146 * t1161 * t1304 + t1147 * t1162 * t1298) * t1311 + ((-t1099 * t1252 - t1101 * t1248 - t1103 * t1244) * MDP(5) + ((t1099 * t1273 + t1101 * t1267 + t1103 * t1261) * MDP(5) + (-t1184 * t1229 - t1185 * t1220 - t1186 * t1211) * MDP(6) + (t1181 * t1229 + t1182 * t1220 + t1183 * t1211) * MDP(7)) * t1187) * t1188; t1098 + (-t1112 * t1209 - t1114 * t1207 - t1116 * t1205) * MDP(6) + (t1112 * t1210 + t1114 * t1208 + t1116 * t1206) * MDP(7) + ((-t1111 * t1252 - t1113 * t1248 - t1115 * t1244) * MDP(5) + ((t1111 * t1273 + t1113 * t1267 + t1115 * t1261) * MDP(5) + (-t1184 * t1230 - t1185 * t1221 - t1186 * t1212) * MDP(6) + (t1181 * t1230 + t1182 * t1221 + t1183 * t1212) * MDP(7)) * t1187) * t1188; (-t1105 * t1253 - t1107 * t1249 - t1109 * t1245) * MDP(6) + (t1105 * t1254 + t1107 * t1250 + t1109 * t1246) * MDP(7) + (-t1151 * t1237 - t1152 * t1228 - t1153 * t1219) * t1311 + ((-t1106 * t1315 - t1108 * t1316 - t1110 * t1317) * MDP(5) + ((t1106 * t1309 + t1108 * t1303 + t1110 * t1297) * MDP(5) + (t1172 * t1231 + t1173 * t1222 + t1174 * t1213) * MDP(6) + (-t1172 * t1232 - t1173 * t1223 - t1174 * t1214) * MDP(7)) * t1187) * t1188; t1098 + (-t1100 * t1253 - t1102 * t1249 - t1104 * t1245) * MDP(6) + (t1100 * t1254 + t1102 * t1250 + t1104 * t1246) * MDP(7) + ((-t1099 * t1315 - t1101 * t1316 - t1103 * t1317) * MDP(5) + ((t1099 * t1309 + t1101 * t1303 + t1103 * t1297) * MDP(5) + (-t1169 * t1231 - t1170 * t1222 - t1171 * t1213) * MDP(6) + (t1169 * t1232 + t1170 * t1223 + t1171 * t1214) * MDP(7)) * t1187) * t1188; (-t1112 * t1253 - t1114 * t1249 - t1116 * t1245) * MDP(6) + (t1112 * t1254 + t1114 * t1250 + t1116 * t1246) * MDP(7) + MDP(8) + (t1148 ^ 2 * t1310 + t1149 ^ 2 * t1304 + t1150 ^ 2 * t1298) * t1311 + ((-t1111 * t1315 - t1113 * t1316 - t1115 * t1317) * MDP(5) + ((t1111 * t1309 + t1113 * t1303 + t1115 * t1297) * MDP(5) + (-t1148 * t1268 - t1149 * t1262 - t1150 * t1256) * MDP(6) + (t1148 * t1269 + t1149 * t1263 + t1150 * t1257) * MDP(7)) * t1187) * t1188;];
%% Postprocessing: Reshape Output
% From vec2mat_3_matlab.m
res = [t1(1), t1(2), t1(3); t1(4), t1(5), t1(6); t1(7), t1(8), t1(9);];
MMX  = res;
