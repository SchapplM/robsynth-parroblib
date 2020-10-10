% Calculate minimal parameter regressor of inertia matrix for parallel robot
% P3PRRRR1G1A0
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
%   pkin=[a2,a4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [12x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3PRRRR1G1A0_convert_par2_MPV_fixb.m

% Output:
% MMX [3x3]
%   minimal parameter regressor of inertia matrix for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 20:34
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMX = P3PRRRR1G1A0_inertia_para_pf_mdp(xP, qJ, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(2,1),zeros(12,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR1G1A0_inertia_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR1G1A0_inertia_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRRR1G1A0_inertia_para_pf_mdp: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR1G1A0_inertia_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR1G1A0_inertia_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [12 1]), ...
  'P3PRRRR1G1A0_inertia_para_pf_mdp: MDP has to be [12x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_MMreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 20:34:16
% EndTime: 2020-03-09 20:34:18
% DurationCPUTime: 1.77s
% Computational Cost: add. (667->205), mult. (2177->425), div. (840->17), fcn. (2742->18), ass. (0->210)
t1243 = legFrame(3,3);
t1216 = sin(t1243);
t1219 = cos(t1243);
t1252 = cos(qJ(3,3));
t1246 = sin(qJ(3,3));
t1253 = cos(qJ(2,3));
t1343 = t1246 * t1253;
t1198 = -t1216 * t1343 - t1219 * t1252;
t1340 = t1252 * t1253;
t1373 = t1216 * t1246;
t1201 = t1219 * t1340 + t1373;
t1395 = t1198 * t1201;
t1202 = -t1216 * t1252 + t1219 * t1343;
t1394 = t1198 * t1202;
t1393 = t1198 * t1216;
t1244 = legFrame(2,3);
t1217 = sin(t1244);
t1220 = cos(t1244);
t1254 = cos(qJ(3,2));
t1248 = sin(qJ(3,2));
t1255 = cos(qJ(2,2));
t1342 = t1248 * t1255;
t1199 = -t1217 * t1342 - t1220 * t1254;
t1339 = t1254 * t1255;
t1372 = t1217 * t1248;
t1203 = t1220 * t1339 + t1372;
t1392 = t1199 * t1203;
t1204 = -t1217 * t1254 + t1220 * t1342;
t1391 = t1199 * t1204;
t1390 = t1199 * t1217;
t1245 = legFrame(1,3);
t1218 = sin(t1245);
t1221 = cos(t1245);
t1256 = cos(qJ(3,1));
t1250 = sin(qJ(3,1));
t1257 = cos(qJ(2,1));
t1341 = t1250 * t1257;
t1200 = -t1218 * t1341 - t1221 * t1256;
t1338 = t1256 * t1257;
t1371 = t1218 * t1250;
t1205 = t1221 * t1338 + t1371;
t1389 = t1200 * t1205;
t1206 = -t1218 * t1256 + t1221 * t1341;
t1388 = t1200 * t1206;
t1387 = t1200 * t1218;
t1370 = t1219 * t1246;
t1207 = t1216 * t1340 - t1370;
t1386 = t1202 * t1207;
t1385 = t1202 * t1219;
t1369 = t1220 * t1248;
t1208 = t1217 * t1339 - t1369;
t1384 = t1204 * t1208;
t1383 = t1204 * t1220;
t1368 = t1221 * t1250;
t1209 = t1218 * t1338 - t1368;
t1382 = t1206 * t1209;
t1381 = t1206 * t1221;
t1231 = 0.1e1 / t1252;
t1235 = 0.1e1 / t1254;
t1239 = 0.1e1 / t1256;
t1260 = t1252 ^ 2;
t1232 = 0.1e1 / t1260;
t1263 = t1254 ^ 2;
t1236 = 0.1e1 / t1263;
t1266 = t1256 ^ 2;
t1240 = 0.1e1 / t1266;
t1259 = 1 / (pkin(2) ^ 2);
t1380 = MDP(2) * t1259;
t1258 = 0.1e1 / pkin(2);
t1379 = (MDP(3) * t1258);
t1378 = (MDP(4) * t1258);
t1377 = MDP(5) * t1259;
t1376 = (MDP(7) * t1259);
t1375 = (MDP(8) * t1259);
t1374 = MDP(9) * t1259;
t1247 = sin(qJ(2,3));
t1223 = 0.1e1 / t1247;
t1367 = t1223 * t1231;
t1366 = t1223 * t1232;
t1233 = t1231 * t1232;
t1365 = t1223 * t1233;
t1364 = t1223 * t1253;
t1224 = 0.1e1 / t1247 ^ 2;
t1363 = t1224 * t1232;
t1362 = t1224 / t1260 ^ 2;
t1361 = t1224 * t1253;
t1249 = sin(qJ(2,2));
t1226 = 0.1e1 / t1249;
t1360 = t1226 * t1235;
t1359 = t1226 * t1236;
t1237 = t1235 * t1236;
t1358 = t1226 * t1237;
t1357 = t1226 * t1255;
t1227 = 0.1e1 / t1249 ^ 2;
t1356 = t1227 * t1236;
t1355 = t1227 / t1263 ^ 2;
t1354 = t1227 * t1255;
t1251 = sin(qJ(2,1));
t1229 = 0.1e1 / t1251;
t1353 = t1229 * t1239;
t1352 = t1229 * t1240;
t1241 = t1239 * t1240;
t1351 = t1229 * t1241;
t1350 = t1229 * t1257;
t1230 = 0.1e1 / t1251 ^ 2;
t1349 = t1230 * t1240;
t1348 = t1230 / t1266 ^ 2;
t1347 = t1230 * t1257;
t1346 = t1233 * t1246;
t1345 = t1237 * t1248;
t1344 = t1241 * t1250;
t1311 = t1231 * t1246 * t1247;
t1296 = t1258 * t1311;
t1305 = t1231 * t1258 * t1364;
t1174 = t1202 * t1305 + t1219 * t1296;
t1310 = t1235 * t1248 * t1249;
t1295 = t1258 * t1310;
t1302 = t1235 * t1258 * t1357;
t1176 = t1204 * t1302 + t1220 * t1295;
t1309 = t1239 * t1250 * t1251;
t1294 = t1258 * t1309;
t1299 = t1239 * t1258 * t1350;
t1178 = t1206 * t1299 + t1221 * t1294;
t1337 = 2 * t1379;
t1336 = 2 * t1378;
t1335 = 2 * MDP(6) * t1259;
t1334 = 2 * t1376;
t1333 = 2 * t1375;
t1332 = t1201 * t1367;
t1331 = t1203 * t1360;
t1330 = t1205 * t1353;
t1329 = t1207 * t1367;
t1328 = t1208 * t1360;
t1327 = t1209 * t1353;
t1222 = t1246 ^ 2;
t1326 = t1222 * t1362;
t1325 = t1232 * t1364;
t1324 = t1223 * t1346;
t1323 = t1224 * t1346;
t1322 = t1233 * t1361;
t1225 = t1248 ^ 2;
t1321 = t1225 * t1355;
t1320 = t1236 * t1357;
t1319 = t1226 * t1345;
t1318 = t1227 * t1345;
t1317 = t1237 * t1354;
t1228 = t1250 ^ 2;
t1316 = t1228 * t1348;
t1315 = t1240 * t1350;
t1314 = t1229 * t1344;
t1313 = t1230 * t1344;
t1312 = t1241 * t1347;
t1308 = t1362 * t1394;
t1307 = t1355 * t1391;
t1306 = t1348 * t1388;
t1304 = t1246 * t1325;
t1303 = t1246 * t1322;
t1301 = t1248 * t1320;
t1300 = t1248 * t1317;
t1298 = t1250 * t1315;
t1297 = t1250 * t1312;
t1293 = t1202 * t1304;
t1292 = t1204 * t1301;
t1291 = t1206 * t1298;
t1290 = t1198 * t1207 + t1201 * t1202;
t1289 = t1199 * t1208 + t1203 * t1204;
t1288 = t1200 * t1209 + t1205 * t1206;
t1287 = t1223 * (-t1198 * t1219 + t1202 * t1216);
t1286 = t1226 * (-t1199 * t1220 + t1204 * t1217);
t1285 = t1229 * (-t1200 * t1221 + t1206 * t1218);
t1284 = t1232 * (t1198 * t1361 - t1373);
t1283 = t1236 * (t1204 * t1354 + t1369);
t1282 = t1240 * (t1200 * t1347 - t1371);
t1281 = t1240 * (t1206 * t1347 + t1368);
t1280 = (t1199 * t1354 - t1372) * t1236;
t1279 = (t1202 * t1361 + t1370) * t1232;
t1278 = -t1198 * t1304 - t1216 * t1247;
t1277 = -t1198 * t1303 - t1216 * t1231;
t1276 = -t1199 * t1301 - t1217 * t1249;
t1275 = -t1199 * t1300 - t1217 * t1235;
t1274 = -t1200 * t1298 - t1218 * t1251;
t1273 = -t1200 * t1297 - t1218 * t1239;
t1272 = -t1202 * t1303 + t1219 * t1231;
t1271 = -t1204 * t1300 + t1220 * t1235;
t1270 = -t1206 * t1297 + t1221 * t1239;
t1269 = (t1288 * t1312 + t1289 * t1317 + t1290 * t1322) * t1379 + (-t1288 * t1351 - t1289 * t1358 - t1290 * t1365) * t1378 + (t1285 * t1344 + t1286 * t1345 + t1287 * t1346) * t1376 + (t1232 * t1287 + t1236 * t1286 + t1240 * t1285) * t1375 + (t1313 * t1388 + t1318 * t1391 + t1323 * t1394) * t1335 + (t1222 * t1308 + t1225 * t1307 + t1228 * t1306) * t1377 + (t1306 + t1307 + t1308) * t1380 + (t1201 * t1207 * t1363 + t1203 * t1208 * t1356 + t1205 * t1209 * t1349) * MDP(1) + (-t1216 * t1219 * t1232 - t1217 * t1220 * t1236 - t1218 * t1221 * t1240) * t1374;
t1215 = t1251 * t1258 * t1221;
t1214 = t1249 * t1258 * t1220;
t1213 = t1247 * t1258 * t1219;
t1197 = t1206 ^ 2;
t1196 = t1204 ^ 2;
t1195 = t1202 ^ 2;
t1194 = t1200 ^ 2;
t1193 = t1199 ^ 2;
t1192 = t1198 ^ 2;
t1190 = t1200 * t1299;
t1188 = t1199 * t1302;
t1186 = t1198 * t1305;
t1185 = -t1258 * t1291 + t1215;
t1184 = t1274 * t1258;
t1183 = -t1258 * t1292 + t1214;
t1182 = t1276 * t1258;
t1181 = -t1258 * t1293 + t1213;
t1180 = t1278 * t1258;
t1177 = -t1218 * t1294 + t1190;
t1175 = -t1217 * t1295 + t1188;
t1173 = -t1216 * t1296 + t1186;
t1164 = (t1327 + t1328 + t1329) * MDP(1) + (t1178 + t1176 + t1174) * MDP(10) + (t1213 + t1214 + t1215) * MDP(11) + ((t1202 * t1325 + t1204 * t1320 + t1206 * t1315) * MDP(3) + (-t1202 * t1232 - t1204 * t1236 - t1206 * t1240) * MDP(4) + (-t1291 - t1292 - t1293) * MDP(11)) * t1258;
t1163 = (t1330 + t1331 + t1332) * MDP(1) + (t1186 + t1188 + t1190) * MDP(10) + ((t1198 * t1325 + t1199 * t1320 + t1200 * t1315) * MDP(3) + (-t1198 * t1232 - t1199 * t1236 - t1200 * t1240) * MDP(4) + (-t1216 * t1311 - t1217 * t1310 - t1218 * t1309) * MDP(10) + (t1274 + t1276 + t1278) * MDP(11)) * t1258;
t1 = [(t1201 ^ 2 * t1363 + t1203 ^ 2 * t1356 + t1205 ^ 2 * t1349) * MDP(1) + (t1192 * t1362 + t1193 * t1355 + t1194 * t1348) * t1380 + (t1312 * t1389 + t1317 * t1392 + t1322 * t1395) * t1337 + (-t1351 * t1389 - t1358 * t1392 - t1365 * t1395) * t1336 + (t1192 * t1326 + t1193 * t1321 + t1194 * t1316) * t1377 + (t1192 * t1323 + t1193 * t1318 + t1194 * t1313) * t1335 + (t1314 * t1387 + t1319 * t1390 + t1324 * t1393) * t1334 + (t1352 * t1387 + t1359 * t1390 + t1366 * t1393) * t1333 + (t1216 ^ 2 * t1232 + t1217 ^ 2 * t1236 + t1218 ^ 2 * t1240) * t1374 + (t1173 * t1332 + t1175 * t1331 + t1177 * t1330 + (t1201 * t1284 + t1203 * t1280 + t1205 * t1282) * t1258) * MDP(10) + (t1180 * t1332 + t1182 * t1331 + t1184 * t1330 + (t1277 * t1201 + t1275 * t1203 + t1273 * t1205) * t1258) * MDP(11) + MDP(12); (t1174 * t1332 + t1176 * t1331 + t1178 * t1330 + (t1207 * t1284 + t1208 * t1280 + t1209 * t1282) * t1258) * MDP(10) + (t1181 * t1332 + t1183 * t1331 + t1185 * t1330 + (t1277 * t1207 + t1275 * t1208 + t1273 * t1209) * t1258) * MDP(11) + t1269; t1163; (t1173 * t1329 + t1175 * t1328 + t1177 * t1327 + (t1201 * t1279 + t1203 * t1283 + t1205 * t1281) * t1258) * MDP(10) + (t1180 * t1329 + t1182 * t1328 + t1184 * t1327 + (t1272 * t1201 + t1271 * t1203 + t1270 * t1205) * t1258) * MDP(11) + t1269; (t1207 ^ 2 * t1363 + t1208 ^ 2 * t1356 + t1209 ^ 2 * t1349) * MDP(1) + (t1195 * t1362 + t1196 * t1355 + t1197 * t1348) * t1380 + (t1312 * t1382 + t1317 * t1384 + t1322 * t1386) * t1337 + (-t1351 * t1382 - t1358 * t1384 - t1365 * t1386) * t1336 + (t1195 * t1326 + t1196 * t1321 + t1197 * t1316) * t1377 + (t1195 * t1323 + t1196 * t1318 + t1197 * t1313) * t1335 + (-t1314 * t1381 - t1319 * t1383 - t1324 * t1385) * t1334 + (-t1352 * t1381 - t1359 * t1383 - t1366 * t1385) * t1333 + (t1219 ^ 2 * t1232 + t1220 ^ 2 * t1236 + t1221 ^ 2 * t1240) * t1374 + (t1174 * t1329 + t1176 * t1328 + t1178 * t1327 + (t1207 * t1279 + t1208 * t1283 + t1209 * t1281) * t1258) * MDP(10) + (t1181 * t1329 + t1183 * t1328 + t1185 * t1327 + (t1272 * t1207 + t1271 * t1208 + t1270 * t1209) * t1258) * MDP(11) + MDP(12); t1164; t1163; t1164; 0.3e1 * MDP(1) + MDP(12);];
%% Postprocessing: Reshape Output
% From vec2mat_3_matlab.m
res = [t1(1), t1(2), t1(3); t1(4), t1(5), t1(6); t1(7), t1(8), t1(9);];
MMX  = res;
