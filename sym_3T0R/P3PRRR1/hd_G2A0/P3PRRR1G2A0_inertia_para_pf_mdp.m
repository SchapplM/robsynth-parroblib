% Calculate minimal parameter regressor of inertia matrix for parallel robot
% P3PRRR1G2A0
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
%   see P3PRRR1G2A0_convert_par2_MPV_fixb.m

% Output:
% MMX [3x3]
%   minimal parameter regressor of inertia matrix for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:18
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMX = P3PRRR1G2A0_inertia_para_pf_mdp(xP, qJ, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(8,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR1G2A0_inertia_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR1G2A0_inertia_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3PRRR1G2A0_inertia_para_pf_mdp: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR1G2A0_inertia_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR1G2A0_inertia_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [8 1]), ...
  'P3PRRR1G2A0_inertia_para_pf_mdp: MDP has to be [8x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_MMreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:18:33
% EndTime: 2020-03-09 21:18:35
% DurationCPUTime: 2.15s
% Computational Cost: add. (4710->214), mult. (2059->423), div. (981->9), fcn. (2856->24), ass. (0->201)
t1194 = legFrame(1,2);
t1185 = sin(t1194);
t1188 = cos(t1194);
t1290 = t1185 * t1188;
t1191 = pkin(7) + qJ(2,1);
t1182 = qJ(3,1) + t1191;
t1167 = sin(t1182);
t1173 = sin(t1191);
t1155 = pkin(2) * t1173 + pkin(3) * t1167;
t1299 = t1155 * t1167;
t1170 = cos(t1182);
t1324 = cos(t1191);
t1149 = -t1324 * t1167 + t1170 * t1173;
t1327 = 0.1e1 / t1149 ^ 2;
t1339 = t1290 * t1299 * t1327;
t1193 = legFrame(2,2);
t1184 = sin(t1193);
t1187 = cos(t1193);
t1291 = t1184 * t1187;
t1190 = pkin(7) + qJ(2,2);
t1181 = qJ(3,2) + t1190;
t1166 = sin(t1181);
t1172 = sin(t1190);
t1154 = pkin(2) * t1172 + pkin(3) * t1166;
t1303 = t1154 * t1166;
t1169 = cos(t1181);
t1325 = cos(t1190);
t1148 = -t1325 * t1166 + t1169 * t1172;
t1328 = 0.1e1 / t1148 ^ 2;
t1338 = t1291 * t1303 * t1328;
t1192 = legFrame(3,2);
t1183 = sin(t1192);
t1186 = cos(t1192);
t1292 = t1183 * t1186;
t1189 = pkin(7) + qJ(2,3);
t1180 = qJ(3,3) + t1189;
t1165 = sin(t1180);
t1171 = sin(t1189);
t1153 = pkin(2) * t1171 + pkin(3) * t1165;
t1307 = t1153 * t1165;
t1168 = cos(t1180);
t1326 = cos(t1189);
t1147 = -t1326 * t1165 + t1168 * t1171;
t1329 = 0.1e1 / t1147 ^ 2;
t1337 = t1292 * t1307 * t1329;
t1162 = t1165 ^ 2;
t1336 = t1162 * t1329;
t1163 = t1166 ^ 2;
t1335 = t1163 * t1328;
t1164 = t1167 ^ 2;
t1334 = t1164 * t1327;
t1195 = sin(qJ(3,3));
t1196 = sin(qJ(3,2));
t1197 = sin(qJ(3,1));
t1198 = cos(qJ(3,3));
t1199 = cos(qJ(3,2));
t1200 = cos(qJ(3,1));
t1333 = (t1198 * t1337 + t1199 * t1338 + t1200 * t1339) * MDP(6) + (-t1195 * t1337 - t1196 * t1338 - t1197 * t1339) * MDP(7);
t1332 = 0.1e1 / t1147;
t1331 = 0.1e1 / t1148;
t1330 = 0.1e1 / t1149;
t1323 = MDP(2) / pkin(2) ^ 2;
t1322 = t1332 ^ 2;
t1321 = t1332 * t1183;
t1320 = t1331 ^ 2;
t1319 = t1331 * t1184;
t1318 = t1330 ^ 2;
t1317 = t1330 * t1185;
t1156 = pkin(2) * t1326 + pkin(3) * t1168;
t1316 = t1332 * t1156;
t1315 = t1332 * t1168;
t1314 = t1329 * t1168;
t1157 = pkin(2) * t1325 + pkin(3) * t1169;
t1313 = t1331 * t1157;
t1312 = t1331 * t1169;
t1311 = t1328 * t1169;
t1158 = pkin(2) * t1324 + pkin(3) * t1170;
t1310 = t1330 * t1158;
t1309 = t1330 * t1170;
t1308 = t1327 * t1170;
t1306 = t1153 * t1183;
t1305 = t1153 * t1186;
t1201 = 0.1e1 / pkin(3);
t1304 = t1153 * t1201;
t1302 = t1154 * t1184;
t1301 = t1154 * t1187;
t1300 = t1154 * t1201;
t1298 = t1155 * t1185;
t1297 = t1155 * t1188;
t1296 = t1155 * t1201;
t1295 = t1165 * t1186;
t1294 = t1166 * t1187;
t1293 = t1167 * t1188;
t1202 = 0.1e1 / pkin(2);
t1289 = t1201 * t1202;
t1288 = (-t1162 * t1292 * t1322 - t1163 * t1291 * t1320 - t1164 * t1290 * t1318) * t1323 + (t1290 + t1291 + t1292) * MDP(1);
t1287 = t1332 * t1305;
t1286 = t1165 * t1321;
t1285 = t1202 * t1321;
t1284 = t1332 * t1289;
t1283 = t1331 * t1301;
t1282 = t1166 * t1319;
t1281 = t1202 * t1319;
t1280 = t1331 * t1289;
t1279 = t1330 * t1297;
t1278 = t1167 * t1317;
t1277 = t1202 * t1317;
t1276 = t1330 * t1289;
t1275 = t1332 * t1306;
t1274 = t1332 * t1295;
t1273 = t1195 * t1315;
t1272 = t1198 * t1315;
t1271 = t1202 * t1315;
t1270 = t1156 * t1314;
t1269 = t1329 * t1295;
t1268 = t1331 * t1302;
t1267 = t1331 * t1294;
t1266 = t1196 * t1312;
t1265 = t1199 * t1312;
t1264 = t1202 * t1312;
t1263 = t1157 * t1311;
t1262 = t1328 * t1294;
t1261 = t1330 * t1298;
t1260 = t1330 * t1293;
t1259 = t1197 * t1309;
t1258 = t1200 * t1309;
t1257 = t1202 * t1309;
t1256 = t1158 * t1308;
t1255 = t1327 * t1293;
t1254 = t1307 * t1322;
t1253 = t1332 * t1286;
t1252 = t1195 * t1286;
t1251 = t1198 * t1286;
t1250 = t1303 * t1320;
t1249 = t1331 * t1282;
t1248 = t1196 * t1282;
t1247 = t1199 * t1282;
t1246 = t1299 * t1318;
t1245 = t1330 * t1278;
t1244 = t1197 * t1278;
t1243 = t1200 * t1278;
t1242 = t1195 * t1274;
t1241 = t1198 * t1274;
t1240 = t1202 * t1274;
t1239 = t1306 * t1314;
t1238 = t1156 * t1269;
t1237 = t1196 * t1267;
t1236 = t1199 * t1267;
t1235 = t1202 * t1267;
t1234 = t1302 * t1311;
t1233 = t1157 * t1262;
t1232 = t1197 * t1260;
t1231 = t1200 * t1260;
t1230 = t1202 * t1260;
t1229 = t1298 * t1308;
t1228 = t1158 * t1255;
t1224 = t1195 * t1254;
t1223 = t1198 * t1254;
t1222 = t1287 * t1315;
t1221 = t1156 * t1253;
t1220 = t1196 * t1250;
t1219 = t1199 * t1250;
t1218 = t1283 * t1312;
t1217 = t1157 * t1249;
t1216 = t1197 * t1246;
t1215 = t1200 * t1246;
t1214 = t1279 * t1309;
t1213 = t1158 * t1245;
t1179 = t1188 ^ 2;
t1178 = t1187 ^ 2;
t1177 = t1186 ^ 2;
t1176 = t1185 ^ 2;
t1175 = t1184 ^ 2;
t1174 = t1183 ^ 2;
t1133 = t1158 * t1276;
t1132 = t1157 * t1280;
t1131 = t1156 * t1284;
t1130 = t1276 * t1297;
t1129 = t1280 * t1301;
t1128 = t1284 * t1305;
t1127 = t1133 - t1257;
t1126 = t1133 - 0.2e1 * t1257;
t1125 = t1132 - t1264;
t1124 = t1132 - 0.2e1 * t1264;
t1123 = t1131 - t1271;
t1122 = t1131 - 0.2e1 * t1271;
t1121 = (t1167 - t1296) * t1277;
t1120 = (0.2e1 * t1167 - t1296) * t1277;
t1119 = (t1166 - t1300) * t1281;
t1118 = (0.2e1 * t1166 - t1300) * t1281;
t1117 = (t1165 - t1304) * t1285;
t1116 = (0.2e1 * t1165 - t1304) * t1285;
t1115 = t1130 - t1230;
t1114 = t1130 - 0.2e1 * t1230;
t1113 = t1129 - t1235;
t1112 = t1129 - 0.2e1 * t1235;
t1111 = t1128 - t1240;
t1110 = t1128 - 0.2e1 * t1240;
t1109 = (t1168 * t1269 + t1169 * t1262 + t1170 * t1255) * t1323;
t1107 = (-t1168 * t1253 - t1169 * t1249 - t1170 * t1245) * t1323;
t1 = [(t1176 + t1175 + t1174) * MDP(1) + (-t1110 * t1241 - t1112 * t1236 - t1114 * t1231) * MDP(6) + (t1110 * t1242 + t1112 * t1237 + t1114 * t1232) * MDP(7) + MDP(8) + (t1177 * t1336 + t1178 * t1335 + t1179 * t1334) * t1323 + ((-t1111 * t1274 - t1113 * t1267 - t1115 * t1260) * MDP(5) + ((t1111 * t1287 + t1113 * t1283 + t1115 * t1279) * MDP(5) + (-t1177 * t1223 - t1178 * t1219 - t1179 * t1215) * MDP(6) + (t1177 * t1224 + t1178 * t1220 + t1179 * t1216) * MDP(7)) * t1201) * t1202; (-t1116 * t1241 - t1118 * t1236 - t1120 * t1231) * MDP(6) + (t1116 * t1242 + t1118 * t1237 + t1120 * t1232) * MDP(7) + ((-t1117 * t1274 - t1119 * t1267 - t1121 * t1260) * MDP(5) + ((t1117 * t1287 + t1119 * t1283 + t1121 * t1279) * MDP(5) + t1333) * t1201) * t1202 + t1288; t1109 + (-t1122 * t1241 - t1124 * t1236 - t1126 * t1231) * MDP(6) + (t1122 * t1242 + t1124 * t1237 + t1126 * t1232) * MDP(7) + ((-t1123 * t1274 - t1125 * t1267 - t1127 * t1260) * MDP(5) + ((t1123 * t1287 + t1125 * t1283 + t1127 * t1279) * MDP(5) + (-t1198 * t1222 - t1199 * t1218 - t1200 * t1214) * MDP(6) + (t1195 * t1222 + t1196 * t1218 + t1197 * t1214) * MDP(7)) * t1201) * t1202; (t1110 * t1251 + t1112 * t1247 + t1114 * t1243) * MDP(6) + (-t1110 * t1252 - t1112 * t1248 - t1114 * t1244) * MDP(7) + ((t1111 * t1286 + t1113 * t1282 + t1115 * t1278) * MDP(5) + ((-t1111 * t1275 - t1113 * t1268 - t1115 * t1261) * MDP(5) + t1333) * t1201) * t1202 + t1288; (t1179 + t1178 + t1177) * MDP(1) + (t1116 * t1251 + t1118 * t1247 + t1120 * t1243) * MDP(6) + (-t1116 * t1252 - t1118 * t1248 - t1120 * t1244) * MDP(7) + MDP(8) + (t1174 * t1336 + t1175 * t1335 + t1176 * t1334) * t1323 + ((t1117 * t1286 + t1119 * t1282 + t1121 * t1278) * MDP(5) + ((-t1117 * t1275 - t1119 * t1268 - t1121 * t1261) * MDP(5) + (-t1174 * t1223 - t1175 * t1219 - t1176 * t1215) * MDP(6) + (t1174 * t1224 + t1175 * t1220 + t1176 * t1216) * MDP(7)) * t1201) * t1202; t1107 + (t1122 * t1251 + t1124 * t1247 + t1126 * t1243) * MDP(6) + (-t1122 * t1252 - t1124 * t1248 - t1126 * t1244) * MDP(7) + ((t1123 * t1286 + t1125 * t1282 + t1127 * t1278) * MDP(5) + ((-t1123 * t1275 - t1125 * t1268 - t1127 * t1261) * MDP(5) + (t1198 * t1239 + t1199 * t1234 + t1200 * t1229) * MDP(6) + (-t1195 * t1239 - t1196 * t1234 - t1197 * t1229) * MDP(7)) * t1201) * t1202; t1109 + (-t1110 * t1272 - t1112 * t1265 - t1114 * t1258) * MDP(6) + (t1110 * t1273 + t1112 * t1266 + t1114 * t1259) * MDP(7) + ((-t1111 * t1315 - t1113 * t1312 - t1115 * t1309) * MDP(5) + ((t1111 * t1316 + t1113 * t1313 + t1115 * t1310) * MDP(5) + (-t1198 * t1238 - t1199 * t1233 - t1200 * t1228) * MDP(6) + (t1195 * t1238 + t1196 * t1233 + t1197 * t1228) * MDP(7)) * t1201) * t1202; t1107 + (-t1116 * t1272 - t1118 * t1265 - t1120 * t1258) * MDP(6) + (t1116 * t1273 + t1118 * t1266 + t1120 * t1259) * MDP(7) + ((-t1117 * t1315 - t1119 * t1312 - t1121 * t1309) * MDP(5) + ((t1117 * t1316 + t1119 * t1313 + t1121 * t1310) * MDP(5) + (t1198 * t1221 + t1199 * t1217 + t1200 * t1213) * MDP(6) + (-t1195 * t1221 - t1196 * t1217 - t1197 * t1213) * MDP(7)) * t1201) * t1202; (-t1122 * t1272 - t1124 * t1265 - t1126 * t1258) * MDP(6) + (t1122 * t1273 + t1124 * t1266 + t1126 * t1259) * MDP(7) + MDP(8) + (t1168 ^ 2 * t1329 + t1169 ^ 2 * t1328 + t1170 ^ 2 * t1327) * t1323 + ((-t1123 * t1315 - t1125 * t1312 - t1127 * t1309) * MDP(5) + ((t1123 * t1316 + t1125 * t1313 + t1127 * t1310) * MDP(5) + (-t1198 * t1270 - t1199 * t1263 - t1200 * t1256) * MDP(6) + (t1195 * t1270 + t1196 * t1263 + t1197 * t1256) * MDP(7)) * t1201) * t1202;];
%% Postprocessing: Reshape Output
% From vec2mat_3_matlab.m
res = [t1(1), t1(2), t1(3); t1(4), t1(5), t1(6); t1(7), t1(8), t1(9);];
MMX  = res;
