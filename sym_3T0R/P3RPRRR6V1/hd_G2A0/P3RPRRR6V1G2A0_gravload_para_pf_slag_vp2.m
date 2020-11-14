% Calculate Gravitation load for parallel robot
% P3RPRRR6V1G2A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2020-08-06 18:37
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RPRRR6V1G2A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR6V1G2A0_gravload_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR6V1G2A0_gravload_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR6V1G2A0_gravload_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR6V1G2A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR6V1G2A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RPRRR6V1G2A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR6V1G2A0_gravload_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR6V1G2A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:36:27
% EndTime: 2020-08-06 18:36:28
% DurationCPUTime: 1.38s
% Computational Cost: add. (780->200), mult. (1089->295), div. (51->10), fcn. (657->59), ass. (0->137)
t1285 = sin(pkin(7));
t1306 = -pkin(6) - pkin(5);
t1232 = t1306 * t1285 - pkin(1);
t1291 = sin(qJ(1,3));
t1218 = t1232 * t1291;
t1286 = cos(pkin(7));
t1297 = cos(qJ(1,3));
t1328 = t1297 * t1306;
t1329 = t1297 * t1285;
t1362 = t1218 - pkin(2) * t1329 - (t1291 * pkin(2) + t1328) * t1286;
t1293 = sin(qJ(1,2));
t1219 = t1232 * t1293;
t1299 = cos(qJ(1,2));
t1326 = t1299 * t1306;
t1327 = t1299 * t1285;
t1361 = t1219 - pkin(2) * t1327 - (t1293 * pkin(2) + t1326) * t1286;
t1295 = sin(qJ(1,1));
t1220 = t1232 * t1295;
t1301 = cos(qJ(1,1));
t1324 = t1301 * t1306;
t1325 = t1301 * t1285;
t1360 = t1220 - pkin(2) * t1325 - (t1295 * pkin(2) + t1324) * t1286;
t1359 = -0.2e1 * pkin(1);
t1358 = -0.2e1 * pkin(2);
t1357 = 0.2e1 * pkin(2);
t1356 = -2 * m(3) * pkin(5) + 2 * mrSges(2,2) - 2 * mrSges(3,3);
t1272 = m(3) * pkin(2) + mrSges(2,1);
t1355 = -0.2e1 * t1272;
t1354 = 0.2e1 * t1306;
t1307 = m(2) + m(3);
t1296 = cos(qJ(3,3));
t1249 = t1296 * pkin(3) + pkin(2);
t1298 = cos(qJ(3,2));
t1250 = t1298 * pkin(3) + pkin(2);
t1300 = cos(qJ(3,1));
t1251 = t1300 * pkin(3) + pkin(2);
t1287 = legFrame(3,2);
t1266 = sin(t1287);
t1269 = cos(t1287);
t1215 = t1269 * g(1) - t1266 * g(2);
t1353 = mrSges(3,1) * t1215;
t1288 = legFrame(2,2);
t1267 = sin(t1288);
t1270 = cos(t1288);
t1216 = t1270 * g(1) - t1267 * g(2);
t1352 = mrSges(3,1) * t1216;
t1289 = legFrame(1,2);
t1268 = sin(t1289);
t1271 = cos(t1289);
t1217 = t1271 * g(1) - t1268 * g(2);
t1351 = mrSges(3,1) * t1217;
t1350 = mrSges(3,2) * t1215;
t1349 = mrSges(3,2) * t1216;
t1348 = mrSges(3,2) * t1217;
t1230 = g(3) * t1356;
t1235 = t1307 * pkin(1) + mrSges(1,1);
t1231 = t1235 * g(3);
t1234 = 0.2e1 * t1272 * g(3);
t1318 = pkin(7) + qJ(3,3);
t1259 = qJ(1,3) + t1318;
t1242 = cos(t1259);
t1321 = -pkin(7) + qJ(3,3);
t1260 = qJ(1,3) - t1321;
t1243 = cos(t1260);
t1276 = qJ(1,3) + pkin(7);
t1253 = sin(t1276);
t1256 = cos(t1276);
t1302 = mrSges(3,2) * g(3);
t1303 = mrSges(1,2) * g(3);
t1304 = mrSges(3,1) * g(3);
t1200 = (-t1302 - t1353) * t1243 / 0.2e1 + (t1304 - t1350) * sin(t1260) / 0.2e1 + (t1302 - t1353) * t1242 / 0.2e1 + (t1304 + t1350) * sin(t1259) / 0.2e1 + (t1215 * t1355 + t1230) * t1256 / 0.2e1 + (t1215 * t1356 + t1234) * t1253 / 0.2e1 + (-t1215 * t1235 + t1303) * t1297 + (t1215 * mrSges(1,2) + t1231) * t1291;
t1347 = t1200 / 0.2e1;
t1319 = pkin(7) + qJ(3,2);
t1261 = qJ(1,2) + t1319;
t1244 = cos(t1261);
t1322 = -pkin(7) + qJ(3,2);
t1262 = qJ(1,2) - t1322;
t1245 = cos(t1262);
t1277 = qJ(1,2) + pkin(7);
t1254 = sin(t1277);
t1257 = cos(t1277);
t1201 = (-t1302 - t1352) * t1245 / 0.2e1 + (t1304 - t1349) * sin(t1262) / 0.2e1 + (t1302 - t1352) * t1244 / 0.2e1 + (t1304 + t1349) * sin(t1261) / 0.2e1 + (t1216 * t1355 + t1230) * t1257 / 0.2e1 + (t1216 * t1356 + t1234) * t1254 / 0.2e1 + (-t1216 * t1235 + t1303) * t1299 + (t1216 * mrSges(1,2) + t1231) * t1293;
t1346 = t1201 / 0.2e1;
t1320 = pkin(7) + qJ(3,1);
t1263 = qJ(1,1) + t1320;
t1246 = cos(t1263);
t1323 = -pkin(7) + qJ(3,1);
t1264 = qJ(1,1) - t1323;
t1247 = cos(t1264);
t1278 = qJ(1,1) + pkin(7);
t1255 = sin(t1278);
t1258 = cos(t1278);
t1202 = (-t1302 - t1351) * t1247 / 0.2e1 + (t1304 - t1348) * sin(t1264) / 0.2e1 + (t1302 - t1351) * t1246 / 0.2e1 + (t1304 + t1348) * sin(t1263) / 0.2e1 + (t1217 * t1355 + t1230) * t1258 / 0.2e1 + (t1217 * t1356 + t1234) * t1255 / 0.2e1 + (-t1217 * t1235 + t1303) * t1301 + (t1217 * mrSges(1,2) + t1231) * t1295;
t1345 = t1202 / 0.2e1;
t1212 = t1266 * g(1) + t1269 * g(2);
t1290 = sin(qJ(3,3));
t1308 = 0.1e1 / pkin(3);
t1344 = (-t1212 * (t1296 * mrSges(3,1) - t1290 * mrSges(3,2)) + (g(3) * t1256 + t1215 * t1253) * (mrSges(3,1) * t1290 + mrSges(3,2) * t1296)) * t1308;
t1213 = t1267 * g(1) + t1270 * g(2);
t1292 = sin(qJ(3,2));
t1343 = (-t1213 * (t1298 * mrSges(3,1) - t1292 * mrSges(3,2)) + (g(3) * t1257 + t1216 * t1254) * (mrSges(3,1) * t1292 + mrSges(3,2) * t1298)) * t1308;
t1214 = t1268 * g(1) + t1271 * g(2);
t1294 = sin(qJ(3,1));
t1342 = (-t1214 * (t1300 * mrSges(3,1) - t1294 * mrSges(3,2)) + (g(3) * t1258 + t1217 * t1255) * (mrSges(3,1) * t1294 + mrSges(3,2) * t1300)) * t1308;
t1341 = t1212 * t1307;
t1340 = t1213 * t1307;
t1339 = t1214 * t1307;
t1338 = t1249 * t1291;
t1337 = t1250 * t1293;
t1336 = t1251 * t1295;
t1335 = t1290 * t1266;
t1334 = t1290 * t1269;
t1333 = t1292 * t1267;
t1332 = t1292 * t1270;
t1331 = t1294 * t1268;
t1330 = t1294 * t1271;
t1317 = (t1291 * t1286 + t1329) * t1296 ^ 2 * pkin(3);
t1316 = (t1293 * t1286 + t1327) * t1298 ^ 2 * pkin(3);
t1315 = (t1295 * t1286 + t1325) * t1300 ^ 2 * pkin(3);
t1314 = ((t1328 + t1338) * t1286 - t1218 + t1249 * t1329) * t1344;
t1313 = ((t1326 + t1337) * t1286 - t1219 + t1250 * t1327) * t1343;
t1312 = ((t1324 + t1336) * t1286 - t1220 + t1251 * t1325) * t1342;
t1281 = 0.1e1 / t1294;
t1280 = 0.1e1 / t1292;
t1279 = 0.1e1 / t1290;
t1265 = t1286 * pkin(1);
t1248 = t1265 + pkin(2);
t1241 = -t1289 + t1278;
t1240 = t1289 + t1278;
t1239 = -t1288 + t1277;
t1238 = t1288 + t1277;
t1237 = -t1287 + t1276;
t1236 = t1287 + t1276;
t1223 = 0.1e1 / (t1265 + t1251);
t1222 = 0.1e1 / (t1265 + t1250);
t1221 = 0.1e1 / (t1265 + t1249);
t1 = [-g(1) * m(4) + ((cos(t1241) + cos(t1240)) * t1345 + (-(t1271 * t1315 + (pkin(3) * t1331 - t1360 * t1271) * t1300 + t1248 * t1331) * t1339 - t1271 * t1312) * t1281) * t1223 + ((cos(t1239) + cos(t1238)) * t1346 + (-(t1270 * t1316 + (pkin(3) * t1333 - t1361 * t1270) * t1298 + t1248 * t1333) * t1340 - t1270 * t1313) * t1280) * t1222 + ((cos(t1237) + cos(t1236)) * t1347 + (-(t1269 * t1317 + (pkin(3) * t1335 - t1362 * t1269) * t1296 + t1248 * t1335) * t1341 - t1269 * t1314) * t1279) * t1221; -g(2) * m(4) + ((-sin(t1240) + sin(t1241)) * t1345 + (-(-t1268 * t1315 + (pkin(3) * t1330 + t1360 * t1268) * t1300 + t1248 * t1330) * t1339 + t1268 * t1312) * t1281) * t1223 + ((-sin(t1238) + sin(t1239)) * t1346 + (-(-t1267 * t1316 + (pkin(3) * t1332 + t1361 * t1267) * t1298 + t1248 * t1332) * t1340 + t1267 * t1313) * t1280) * t1222 + ((-sin(t1236) + sin(t1237)) * t1347 + (-(-t1266 * t1317 + (pkin(3) * t1334 + t1362 * t1266) * t1296 + t1248 * t1334) * t1341 + t1266 * t1314) * t1279) * t1221; (t1255 * t1354 + t1301 * t1359 + t1258 * t1358 + (-t1246 - t1247) * pkin(3)) / (pkin(3) * sin(0.2e1 * qJ(3,1)) + t1294 * t1357 + (sin(t1320) + sin(t1323)) * pkin(1)) * t1342 + (t1254 * t1354 + t1299 * t1359 + t1257 * t1358 + (-t1244 - t1245) * pkin(3)) / (pkin(3) * sin(0.2e1 * qJ(3,2)) + t1292 * t1357 + (sin(t1319) + sin(t1322)) * pkin(1)) * t1343 + (t1253 * t1354 + t1297 * t1359 + t1256 * t1358 + (-t1242 - t1243) * pkin(3)) / (pkin(3) * sin(0.2e1 * qJ(3,3)) + t1290 * t1357 + (sin(t1318) + sin(t1321)) * pkin(1)) * t1344 - g(3) * m(4) + (-t1255 * t1202 - ((t1251 * t1301 - t1295 * t1306) * t1286 - t1232 * t1301 - t1285 * t1336) * t1300 * t1281 * t1339) * t1223 + (-t1254 * t1201 - ((t1250 * t1299 - t1293 * t1306) * t1286 - t1232 * t1299 - t1285 * t1337) * t1298 * t1280 * t1340) * t1222 + (-t1253 * t1200 - ((t1249 * t1297 - t1291 * t1306) * t1286 - t1232 * t1297 - t1285 * t1338) * t1296 * t1279 * t1341) * t1221;];
taugX  = t1;
