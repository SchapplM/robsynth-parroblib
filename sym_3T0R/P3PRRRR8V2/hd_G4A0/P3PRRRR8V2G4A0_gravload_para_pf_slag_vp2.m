% Calculate Gravitation load for parallel robot
% P3PRRRR8V2G4A0
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
% Datum: 2020-08-06 18:17
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3PRRRR8V2G4A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(8,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V2G4A0_gravload_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V2G4A0_gravload_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P3PRRRR8V2G4A0_gravload_para_pf_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V2G4A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR8V2G4A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRRR8V2G4A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V2G4A0_gravload_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V2G4A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:15:06
% EndTime: 2020-08-06 18:15:09
% DurationCPUTime: 2.67s
% Computational Cost: add. (1941->263), mult. (4251->510), div. (36->10), fcn. (4494->34), ass. (0->197)
t1294 = sin(qJ(2,3));
t1300 = cos(qJ(2,3));
t1284 = legFrame(3,3);
t1255 = sin(t1284);
t1287 = legFrame(3,1);
t1258 = sin(t1287);
t1261 = cos(t1284);
t1264 = cos(t1287);
t1290 = legFrame(3,2);
t1267 = sin(t1290);
t1346 = t1267 * t1264;
t1349 = t1258 * t1267;
t1270 = cos(t1290);
t1374 = g(1) * t1270;
t1179 = -t1255 * t1374 + (-t1255 * t1349 + t1261 * t1264) * g(2) + (t1255 * t1346 + t1258 * t1261) * g(3);
t1180 = t1261 * t1374 + (t1255 * t1264 + t1261 * t1349) * g(2) + (t1255 * t1258 - t1261 * t1346) * g(3);
t1280 = sin(pkin(8));
t1282 = cos(pkin(8));
t1321 = t1179 * t1280 + t1180 * t1282;
t1212 = g(1) * t1267 + (-g(2) * t1258 + g(3) * t1264) * t1270;
t1281 = sin(pkin(4));
t1283 = cos(pkin(4));
t1322 = t1179 * t1282 - t1180 * t1280;
t1382 = t1212 * t1281 + t1322 * t1283;
t1310 = t1382 * t1294 + t1321 * t1300;
t1296 = sin(qJ(2,2));
t1302 = cos(qJ(2,2));
t1285 = legFrame(2,3);
t1256 = sin(t1285);
t1288 = legFrame(2,1);
t1259 = sin(t1288);
t1262 = cos(t1285);
t1265 = cos(t1288);
t1291 = legFrame(2,2);
t1268 = sin(t1291);
t1345 = t1268 * t1265;
t1348 = t1259 * t1268;
t1271 = cos(t1291);
t1373 = g(1) * t1271;
t1181 = -t1256 * t1373 + (-t1256 * t1348 + t1262 * t1265) * g(2) + (t1256 * t1345 + t1259 * t1262) * g(3);
t1182 = t1262 * t1373 + (t1256 * t1265 + t1262 * t1348) * g(2) + (t1256 * t1259 - t1262 * t1345) * g(3);
t1319 = t1181 * t1280 + t1182 * t1282;
t1213 = g(1) * t1268 + (-g(2) * t1259 + g(3) * t1265) * t1271;
t1320 = t1181 * t1282 - t1182 * t1280;
t1383 = t1213 * t1281 + t1320 * t1283;
t1309 = t1383 * t1296 + t1319 * t1302;
t1298 = sin(qJ(2,1));
t1304 = cos(qJ(2,1));
t1286 = legFrame(1,3);
t1257 = sin(t1286);
t1289 = legFrame(1,1);
t1260 = sin(t1289);
t1263 = cos(t1286);
t1266 = cos(t1289);
t1292 = legFrame(1,2);
t1269 = sin(t1292);
t1344 = t1269 * t1266;
t1347 = t1260 * t1269;
t1272 = cos(t1292);
t1372 = g(1) * t1272;
t1183 = -t1257 * t1372 + (-t1257 * t1347 + t1263 * t1266) * g(2) + (t1257 * t1344 + t1260 * t1263) * g(3);
t1184 = t1263 * t1372 + (t1257 * t1266 + t1263 * t1347) * g(2) + (t1257 * t1260 - t1263 * t1344) * g(3);
t1317 = t1183 * t1280 + t1184 * t1282;
t1214 = g(1) * t1269 + (-g(2) * t1260 + g(3) * t1266) * t1272;
t1318 = t1183 * t1282 - t1184 * t1280;
t1384 = t1214 * t1281 + t1318 * t1283;
t1308 = t1384 * t1298 + t1317 * t1304;
t1221 = t1255 * t1282 + t1261 * t1280;
t1393 = t1221 * t1281;
t1222 = t1256 * t1282 + t1262 * t1280;
t1392 = t1222 * t1281;
t1223 = t1257 * t1282 + t1263 * t1280;
t1391 = t1223 * t1281;
t1387 = -t1214 * t1283 + t1318 * t1281;
t1386 = -t1213 * t1283 + t1320 * t1281;
t1385 = -t1212 * t1283 + t1322 * t1281;
t1299 = cos(qJ(3,3));
t1381 = pkin(3) * t1299 ^ 2;
t1301 = cos(qJ(3,2));
t1380 = pkin(3) * t1301 ^ 2;
t1303 = cos(qJ(3,1));
t1379 = pkin(3) * t1303 ^ 2;
t1378 = pkin(3) * t1281;
t1377 = pkin(3) * t1299;
t1376 = pkin(3) * t1301;
t1375 = pkin(3) * t1303;
t1293 = sin(qJ(3,3));
t1371 = t1293 * pkin(2);
t1295 = sin(qJ(3,2));
t1370 = t1295 * pkin(2);
t1297 = sin(qJ(3,1));
t1369 = t1297 * pkin(2);
t1368 = m(3) * pkin(2) + mrSges(2,1);
t1248 = pkin(2) + t1377;
t1305 = pkin(7) + pkin(6);
t1251 = t1305 * t1300;
t1233 = t1248 * t1294 - t1251;
t1337 = t1283 * t1293;
t1239 = t1248 * t1337;
t1340 = t1281 * t1299;
t1367 = ((t1385 * mrSges(3,1) + t1310 * mrSges(3,2)) * t1299 + (t1310 * mrSges(3,1) - t1385 * mrSges(3,2)) * t1293) / (t1233 * t1340 + t1239);
t1249 = pkin(2) + t1376;
t1252 = t1305 * t1302;
t1234 = t1249 * t1296 - t1252;
t1335 = t1283 * t1295;
t1240 = t1249 * t1335;
t1339 = t1281 * t1301;
t1366 = ((t1386 * mrSges(3,1) + t1309 * mrSges(3,2)) * t1301 + (t1309 * mrSges(3,1) - t1386 * mrSges(3,2)) * t1295) / (t1234 * t1339 + t1240);
t1250 = pkin(2) + t1375;
t1253 = t1305 * t1304;
t1235 = t1250 * t1298 - t1253;
t1333 = t1283 * t1297;
t1241 = t1250 * t1333;
t1338 = t1281 * t1303;
t1365 = ((t1387 * mrSges(3,1) + t1308 * mrSges(3,2)) * t1303 + (t1308 * mrSges(3,1) - t1387 * mrSges(3,2)) * t1297) / (t1235 * t1338 + t1241);
t1254 = m(3) * pkin(6) - mrSges(2,2) + mrSges(3,3);
t1164 = -t1310 * t1254 + (t1321 * t1294 - t1300 * t1382) * (mrSges(3,1) * t1299 - mrSges(3,2) * t1293 + t1368);
t1242 = pkin(2) * t1294 - t1251;
t1364 = t1164 / (t1239 + (t1294 * t1377 + t1242) * t1340);
t1165 = -t1309 * t1254 + (t1319 * t1296 - t1302 * t1383) * (mrSges(3,1) * t1301 - mrSges(3,2) * t1295 + t1368);
t1243 = pkin(2) * t1296 - t1252;
t1363 = t1165 / (t1240 + (t1296 * t1376 + t1243) * t1339);
t1166 = -t1308 * t1254 + (t1317 * t1298 - t1304 * t1384) * (mrSges(3,1) * t1303 - mrSges(3,2) * t1297 + t1368);
t1244 = pkin(2) * t1298 - t1253;
t1362 = t1166 / (t1241 + (t1298 * t1375 + t1244) * t1338);
t1215 = pkin(3) * t1337 + t1242 * t1281;
t1331 = t1294 * t1281;
t1191 = 0.1e1 / (pkin(2) * t1337 + t1215 * t1299 + t1331 * t1381);
t1361 = t1191 * t1212;
t1216 = pkin(3) * t1335 + t1243 * t1281;
t1329 = t1296 * t1281;
t1192 = 0.1e1 / (pkin(2) * t1335 + t1216 * t1301 + t1329 * t1380);
t1360 = t1192 * t1213;
t1217 = pkin(3) * t1333 + t1244 * t1281;
t1327 = t1298 * t1281;
t1193 = 0.1e1 / (pkin(2) * t1333 + t1217 * t1303 + t1327 * t1379);
t1359 = t1193 * t1214;
t1330 = t1294 * t1305;
t1352 = (t1248 * t1300 + t1330) * t1283;
t1328 = t1296 * t1305;
t1351 = (t1249 * t1302 + t1328) * t1283;
t1326 = t1298 * t1305;
t1350 = (t1250 * t1304 + t1326) * t1283;
t1218 = -t1255 * t1280 + t1261 * t1282;
t1343 = t1281 * t1218;
t1219 = -t1256 * t1280 + t1262 * t1282;
t1342 = t1281 * t1219;
t1220 = -t1257 * t1280 + t1263 * t1282;
t1341 = t1281 * t1220;
t1336 = t1283 * t1294;
t1334 = t1283 * t1296;
t1332 = t1283 * t1298;
t1224 = t1280 * t1336 - t1282 * t1300;
t1227 = t1280 * t1300 + t1282 * t1336;
t1316 = t1224 * t1261 + t1227 * t1255;
t1225 = t1280 * t1334 - t1282 * t1302;
t1228 = t1280 * t1302 + t1282 * t1334;
t1315 = t1225 * t1262 + t1228 * t1256;
t1226 = t1280 * t1332 - t1282 * t1304;
t1229 = t1280 * t1304 + t1282 * t1332;
t1314 = t1226 * t1263 + t1229 * t1257;
t1313 = -t1242 * t1283 + t1293 * t1378;
t1312 = -t1243 * t1283 + t1295 * t1378;
t1311 = -t1244 * t1283 + t1297 * t1378;
t1307 = 0.1e1 / pkin(3);
t1276 = m(1) + m(2) + m(3);
t1247 = pkin(2) * t1304 + t1326;
t1246 = pkin(2) * t1302 + t1328;
t1245 = pkin(2) * t1300 + t1330;
t1205 = -t1269 * t1391 + t1272 * t1283;
t1204 = -t1268 * t1392 + t1271 * t1283;
t1203 = -t1267 * t1393 + t1270 * t1283;
t1202 = t1247 * t1280 - t1311 * t1282;
t1201 = t1246 * t1280 - t1312 * t1282;
t1200 = t1245 * t1280 - t1313 * t1282;
t1199 = -t1247 * t1282 - t1311 * t1280;
t1198 = -t1246 * t1282 - t1312 * t1280;
t1197 = -t1245 * t1282 - t1313 * t1280;
t1196 = -t1226 * t1257 + t1229 * t1263;
t1195 = -t1225 * t1256 + t1228 * t1262;
t1194 = -t1224 * t1255 + t1227 * t1261;
t1190 = t1220 * t1347 + t1223 * t1266;
t1189 = t1219 * t1348 + t1222 * t1265;
t1188 = t1218 * t1349 + t1221 * t1264;
t1187 = t1220 * t1344 - t1223 * t1260;
t1186 = t1219 * t1345 - t1222 * t1259;
t1185 = t1218 * t1346 - t1221 * t1258;
t1178 = t1314 * t1269 + t1272 * t1327;
t1177 = t1315 * t1268 + t1271 * t1329;
t1176 = t1316 * t1267 + t1270 * t1331;
t1175 = -t1199 * t1257 + t1202 * t1263;
t1174 = -t1198 * t1256 + t1201 * t1262;
t1173 = -t1197 * t1255 + t1200 * t1261;
t1169 = t1217 * t1272 + (t1199 * t1263 + t1202 * t1257) * t1269;
t1168 = t1216 * t1271 + (t1198 * t1262 + t1201 * t1256) * t1268;
t1167 = t1215 * t1270 + (t1197 * t1261 + t1200 * t1255) * t1267;
t1 = [-t1272 * (t1220 * t1338 + t1297 * (t1220 * t1332 + t1223 * t1304)) * t1193 * t1166 - t1271 * (t1219 * t1339 + t1295 * (t1219 * t1334 + t1222 * t1302)) * t1192 * t1165 - t1270 * (t1218 * t1340 + t1293 * (t1218 * t1336 + t1221 * t1300)) * t1191 * t1164 - g(1) * m(4) + (t1272 * (-t1220 * t1350 + t1223 * t1235) * t1365 + t1271 * (-t1219 * t1351 + t1222 * t1234) * t1366 + t1270 * (-t1218 * t1352 + t1221 * t1233) * t1367) * t1307 + (-(-((-t1220 * t1304 + t1223 * t1332) * t1272 - t1269 * t1327) * t1379 + ((t1220 * t1247 + t1311 * t1223) * t1272 + t1217 * t1269) * t1303 + (t1269 * t1283 + t1272 * t1391) * t1369) * t1359 - (-((-t1219 * t1302 + t1222 * t1334) * t1271 - t1268 * t1329) * t1380 + ((t1219 * t1246 + t1312 * t1222) * t1271 + t1216 * t1268) * t1301 + (t1268 * t1283 + t1271 * t1392) * t1370) * t1360 - (-((-t1218 * t1300 + t1221 * t1336) * t1270 - t1267 * t1331) * t1381 + ((t1218 * t1245 + t1313 * t1221) * t1270 + t1215 * t1267) * t1299 + (t1267 * t1283 + t1270 * t1393) * t1371) * t1361) * t1276; ((-t1196 * t1347 - t1314 * t1266) * t1297 - t1190 * t1338) * t1362 + ((-t1195 * t1348 - t1315 * t1265) * t1295 - t1189 * t1339) * t1363 + ((-t1194 * t1349 - t1316 * t1264) * t1293 - t1188 * t1340) * t1364 - g(2) * m(4) + ((-t1190 * t1350 + (-t1220 * t1266 + t1223 * t1347) * t1235) * t1365 + (-t1189 * t1351 + (-t1219 * t1265 + t1222 * t1348) * t1234) * t1366 + (-t1188 * t1352 + (-t1218 * t1264 + t1221 * t1349) * t1233) * t1367) * t1307 + (-(-(t1178 * t1260 - t1196 * t1266) * t1379 + (-t1169 * t1260 + t1175 * t1266) * t1303 - (t1205 * t1260 + t1266 * t1341) * t1369) * t1359 - (-(t1177 * t1259 - t1195 * t1265) * t1380 + (-t1168 * t1259 + t1174 * t1265) * t1301 - (t1204 * t1259 + t1265 * t1342) * t1370) * t1360 - (-(t1176 * t1258 - t1194 * t1264) * t1381 + (-t1167 * t1258 + t1173 * t1264) * t1299 - (t1203 * t1258 + t1264 * t1343) * t1371) * t1361) * t1276; ((t1196 * t1344 - t1260 * t1314) * t1297 + t1187 * t1338) * t1362 + ((t1195 * t1345 - t1259 * t1315) * t1295 + t1186 * t1339) * t1363 + ((t1194 * t1346 - t1258 * t1316) * t1293 + t1185 * t1340) * t1364 - g(3) * m(4) + ((t1187 * t1350 - (t1220 * t1260 + t1223 * t1344) * t1235) * t1365 + (t1186 * t1351 - (t1219 * t1259 + t1222 * t1345) * t1234) * t1366 + (t1185 * t1352 - (t1218 * t1258 + t1221 * t1346) * t1233) * t1367) * t1307 + (-((t1178 * t1266 + t1196 * t1260) * t1379 + (t1169 * t1266 + t1175 * t1260) * t1303 + (t1205 * t1266 - t1260 * t1341) * t1369) * t1359 - ((t1177 * t1265 + t1195 * t1259) * t1380 + (t1168 * t1265 + t1174 * t1259) * t1301 + (t1204 * t1265 - t1259 * t1342) * t1370) * t1360 - ((t1176 * t1264 + t1194 * t1258) * t1381 + (t1167 * t1264 + t1173 * t1258) * t1299 + (t1203 * t1264 - t1258 * t1343) * t1371) * t1361) * t1276;];
taugX  = t1;
