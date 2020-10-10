% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3PRRRR8V1G1A0
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4,alpha2,d2,d4,theta1]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3x12]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 16:50
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3PRRRR8V1G1A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V1G1A0_gravload_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V1G1A0_gravload_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR8V1G1A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3PRRRR8V1G1A0_gravload_para_pf_regmin: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V1G1A0_gravload_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V1G1A0_gravload_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 16:50:26
% EndTime: 2020-08-06 16:50:28
% DurationCPUTime: 1.15s
% Computational Cost: add. (653->175), mult. (1621->354), div. (84->7), fcn. (1908->22), ass. (0->150)
t1280 = cos(qJ(3,1));
t1281 = cos(qJ(2,1));
t1275 = sin(qJ(2,1));
t1316 = t1275 * t1280;
t1244 = pkin(2) * t1316 - t1281 * pkin(5);
t1264 = sin(pkin(3));
t1266 = cos(pkin(3));
t1274 = sin(qJ(3,1));
t1338 = t1266 * t1274;
t1289 = 0.1e1 / (pkin(2) * t1338 + t1244 * t1264);
t1358 = t1289 / t1280;
t1278 = cos(qJ(3,2));
t1279 = cos(qJ(2,2));
t1273 = sin(qJ(2,2));
t1322 = t1273 * t1278;
t1243 = pkin(2) * t1322 - t1279 * pkin(5);
t1272 = sin(qJ(3,2));
t1340 = t1266 * t1272;
t1290 = 0.1e1 / (pkin(2) * t1340 + t1243 * t1264);
t1357 = t1290 / t1278;
t1276 = cos(qJ(3,3));
t1277 = cos(qJ(2,3));
t1271 = sin(qJ(2,3));
t1328 = t1271 * t1276;
t1242 = pkin(2) * t1328 - t1277 * pkin(5);
t1270 = sin(qJ(3,3));
t1342 = t1266 * t1270;
t1291 = 0.1e1 / (pkin(2) * t1342 + t1242 * t1264);
t1356 = t1291 / t1276;
t1355 = pkin(2) * t1276;
t1354 = pkin(2) * t1278;
t1353 = pkin(2) * t1280;
t1352 = t1264 * g(3);
t1263 = sin(pkin(6));
t1265 = cos(pkin(6));
t1248 = t1263 * g(1) - t1265 * g(2);
t1249 = t1265 * g(1) + t1263 * g(2);
t1267 = legFrame(3,3);
t1254 = sin(t1267);
t1257 = cos(t1267);
t1191 = (-t1352 + (t1248 * t1257 + t1249 * t1254) * t1266) * t1277 + (-t1254 * t1248 + t1249 * t1257) * t1271;
t1351 = t1191 * t1291;
t1268 = legFrame(2,3);
t1255 = sin(t1268);
t1258 = cos(t1268);
t1192 = (-t1352 + (t1248 * t1258 + t1249 * t1255) * t1266) * t1279 + (-t1255 * t1248 + t1249 * t1258) * t1273;
t1350 = t1192 * t1290;
t1269 = legFrame(1,3);
t1256 = sin(t1269);
t1259 = cos(t1269);
t1193 = (-t1352 + (t1248 * t1259 + t1249 * t1256) * t1266) * t1281 + (-t1256 * t1248 + t1249 * t1259) * t1275;
t1349 = t1193 * t1289;
t1341 = t1266 * t1271;
t1339 = t1266 * t1273;
t1337 = t1266 * t1275;
t1336 = t1266 * t1277;
t1335 = t1266 * t1279;
t1334 = t1266 * t1281;
t1333 = t1270 * t1264;
t1332 = t1270 * t1271;
t1331 = t1270 * t1277;
t1218 = -t1263 * t1254 + t1257 * t1265;
t1330 = t1271 * t1218;
t1221 = t1265 * t1254 + t1257 * t1263;
t1329 = t1271 * t1221;
t1327 = t1272 * t1264;
t1326 = t1272 * t1273;
t1325 = t1272 * t1279;
t1219 = -t1263 * t1255 + t1258 * t1265;
t1324 = t1273 * t1219;
t1222 = t1265 * t1255 + t1258 * t1263;
t1323 = t1273 * t1222;
t1321 = t1274 * t1264;
t1320 = t1274 * t1275;
t1319 = t1274 * t1281;
t1220 = -t1263 * t1256 + t1259 * t1265;
t1318 = t1275 * t1220;
t1223 = t1265 * t1256 + t1259 * t1263;
t1317 = t1275 * t1223;
t1315 = t1276 * t1264;
t1314 = t1276 * t1277;
t1313 = t1277 * t1218;
t1312 = t1277 * t1221;
t1311 = t1278 * t1264;
t1310 = t1278 * t1279;
t1309 = t1279 * t1219;
t1308 = t1279 * t1222;
t1307 = t1280 * t1264;
t1306 = t1280 * t1281;
t1305 = t1281 * t1220;
t1304 = t1281 * t1223;
t1303 = (-(t1266 * t1313 - t1329) * t1355 - pkin(5) * (t1266 * t1330 + t1312)) * t1356;
t1302 = (-(t1266 * t1312 + t1330) * t1355 - (t1266 * t1329 - t1313) * pkin(5)) * t1356;
t1301 = (-(t1266 * t1309 - t1323) * t1354 - pkin(5) * (t1266 * t1324 + t1308)) * t1357;
t1300 = (-(t1266 * t1308 + t1324) * t1354 - (t1266 * t1323 - t1309) * pkin(5)) * t1357;
t1299 = (-(t1266 * t1305 - t1317) * t1353 - pkin(5) * (t1266 * t1318 + t1304)) * t1358;
t1298 = (-(t1266 * t1304 + t1318) * t1353 - (t1266 * t1317 - t1305) * pkin(5)) * t1358;
t1236 = t1254 * g(1) - t1257 * g(2);
t1239 = t1257 * g(1) + t1254 * g(2);
t1297 = (t1239 * (t1263 * t1336 + t1265 * t1271) + t1236 * (-t1263 * t1271 + t1265 * t1336) - t1277 * t1352) * t1356;
t1225 = t1263 * t1277 + t1265 * t1341;
t1228 = -t1263 * t1341 + t1265 * t1277;
t1296 = (-t1236 * t1225 + t1239 * t1228 + t1271 * t1352) * t1356;
t1237 = t1255 * g(1) - t1258 * g(2);
t1240 = t1258 * g(1) + t1255 * g(2);
t1295 = (t1240 * (t1263 * t1335 + t1265 * t1273) + t1237 * (-t1263 * t1273 + t1265 * t1335) - t1279 * t1352) * t1357;
t1226 = t1263 * t1279 + t1265 * t1339;
t1229 = -t1263 * t1339 + t1265 * t1279;
t1294 = (-t1237 * t1226 + t1240 * t1229 + t1273 * t1352) * t1357;
t1238 = t1256 * g(1) - t1259 * g(2);
t1241 = t1259 * g(1) + t1256 * g(2);
t1293 = (t1241 * (t1263 * t1334 + t1265 * t1275) + t1238 * (-t1263 * t1275 + t1265 * t1334) - t1281 * t1352) * t1358;
t1224 = t1263 * t1337 - t1265 * t1281;
t1227 = t1263 * t1281 + t1265 * t1337;
t1292 = (-t1241 * t1224 - t1238 * t1227 + t1275 * t1352) * t1358;
t1288 = t1191 * t1270 * t1356;
t1287 = t1192 * t1272 * t1357;
t1286 = t1193 * t1274 * t1358;
t1285 = pkin(2) * t1333 - t1242 * t1266;
t1284 = pkin(2) * t1327 - t1243 * t1266;
t1283 = pkin(2) * t1321 - t1244 * t1266;
t1282 = 0.1e1 / pkin(2);
t1247 = pkin(2) * t1306 + pkin(5) * t1275;
t1246 = pkin(2) * t1310 + pkin(5) * t1273;
t1245 = pkin(2) * t1314 + pkin(5) * t1271;
t1235 = t1266 * t1316 - t1321;
t1234 = t1266 * t1320 + t1307;
t1233 = t1266 * t1322 - t1327;
t1232 = t1266 * t1328 - t1333;
t1231 = t1266 * t1326 + t1311;
t1230 = t1266 * t1332 + t1315;
t1211 = t1263 * t1247 - t1283 * t1265;
t1210 = t1263 * t1246 - t1284 * t1265;
t1209 = t1263 * t1245 - t1285 * t1265;
t1208 = t1247 * t1265 + t1283 * t1263;
t1207 = t1246 * t1265 + t1284 * t1263;
t1206 = t1245 * t1265 + t1285 * t1263;
t1199 = (-t1224 * t1259 - t1256 * t1227) * t1274 - t1223 * t1307;
t1198 = (-t1255 * t1226 + t1229 * t1258) * t1272 - t1222 * t1311;
t1197 = (-t1254 * t1225 + t1228 * t1257) * t1270 - t1221 * t1315;
t1196 = (t1256 * t1224 - t1227 * t1259) * t1274 - t1220 * t1307;
t1195 = (-t1226 * t1258 - t1255 * t1229) * t1272 - t1219 * t1311;
t1194 = (-t1225 * t1257 - t1254 * t1228) * t1270 - t1218 * t1315;
t1184 = (-t1235 * t1263 + t1265 * t1306) * t1241 - t1238 * (t1235 * t1265 + t1263 * t1306) + g(3) * (t1275 * t1307 + t1338);
t1183 = (-t1233 * t1263 + t1265 * t1310) * t1240 - t1237 * (t1233 * t1265 + t1263 * t1310) + g(3) * (t1273 * t1311 + t1340);
t1182 = (-t1232 * t1263 + t1265 * t1314) * t1239 - t1236 * (t1232 * t1265 + t1263 * t1314) + g(3) * (t1271 * t1315 + t1342);
t1181 = t1241 * (-t1234 * t1263 + t1265 * t1319) - (t1234 * t1265 + t1263 * t1319) * t1238 - g(3) * (-t1264 * t1320 + t1266 * t1280);
t1180 = t1240 * (-t1231 * t1263 + t1265 * t1325) - (t1231 * t1265 + t1263 * t1325) * t1237 - g(3) * (-t1264 * t1326 + t1266 * t1278);
t1179 = t1239 * (-t1230 * t1263 + t1265 * t1331) - (t1230 * t1265 + t1263 * t1331) * t1236 - g(3) * (-t1264 * t1332 + t1266 * t1276);
t1 = [(-(t1208 * t1259 - t1256 * t1211) * t1289 - (t1207 * t1258 - t1255 * t1210) * t1290 - (t1206 * t1257 - t1254 * t1209) * t1291) * g(3), 0, t1194 * t1297 + t1195 * t1295 + t1196 * t1293, t1194 * t1296 + t1195 * t1294 + t1196 * t1292, 0, 0, 0, 0, 0, t1194 * t1351 + t1195 * t1350 + t1196 * t1349 + (t1179 * t1303 + t1180 * t1301 + t1181 * t1299) * t1282, -t1194 * t1288 - t1195 * t1287 - t1196 * t1286 + (t1182 * t1303 + t1183 * t1301 + t1184 * t1299) * t1282, -g(1); (-(t1208 * t1256 + t1211 * t1259) * t1289 - (t1207 * t1255 + t1210 * t1258) * t1290 - (t1206 * t1254 + t1209 * t1257) * t1291) * g(3), 0, t1197 * t1297 + t1198 * t1295 + t1199 * t1293, t1197 * t1296 + t1198 * t1294 + t1199 * t1292, 0, 0, 0, 0, 0, t1197 * t1351 + t1198 * t1350 + t1199 * t1349 + (t1179 * t1302 + t1180 * t1300 + t1181 * t1298) * t1282, -t1197 * t1288 - t1198 * t1287 - t1199 * t1286 + (t1182 * t1302 + t1183 * t1300 + t1184 * t1298) * t1282, -g(2); -0.3e1 * g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3);];
tau_reg  = t1;
