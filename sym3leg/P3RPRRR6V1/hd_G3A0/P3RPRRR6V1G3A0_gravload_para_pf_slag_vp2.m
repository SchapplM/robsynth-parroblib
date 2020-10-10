% Calculate Gravitation load for parallel robot
% P3RPRRR6V1G3A0
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
% Datum: 2020-08-06 18:43
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RPRRR6V1G3A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR6V1G3A0_gravload_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR6V1G3A0_gravload_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR6V1G3A0_gravload_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR6V1G3A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR6V1G3A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RPRRR6V1G3A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR6V1G3A0_gravload_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR6V1G3A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:41:53
% EndTime: 2020-08-06 18:41:54
% DurationCPUTime: 1.21s
% Computational Cost: add. (780->200), mult. (1089->293), div. (51->10), fcn. (657->59), ass. (0->136)
t1273 = sin(pkin(7));
t1293 = -pkin(6) - pkin(5);
t1219 = t1293 * t1273 - pkin(1);
t1285 = cos(qJ(1,3));
t1206 = t1219 * t1285;
t1274 = cos(pkin(7));
t1279 = sin(qJ(1,3));
t1319 = t1279 * t1293;
t1320 = t1279 * t1273;
t1346 = t1206 + pkin(2) * t1320 - (pkin(2) * t1285 - t1319) * t1274;
t1287 = cos(qJ(1,2));
t1207 = t1219 * t1287;
t1281 = sin(qJ(1,2));
t1315 = t1281 * t1293;
t1316 = t1281 * t1273;
t1345 = t1207 + pkin(2) * t1316 - (pkin(2) * t1287 - t1315) * t1274;
t1289 = cos(qJ(1,1));
t1208 = t1219 * t1289;
t1283 = sin(qJ(1,1));
t1311 = t1283 * t1293;
t1312 = t1283 * t1273;
t1344 = t1208 + pkin(2) * t1312 - (pkin(2) * t1289 - t1311) * t1274;
t1343 = 0.2e1 * pkin(1);
t1342 = 0.2e1 * pkin(2);
t1240 = (-m(3) * pkin(5) + mrSges(2,2) - mrSges(3,3));
t1341 = 2 * t1240;
t1340 = 0.2e1 * m(3) * pkin(2) + (2 * mrSges(2,1));
t1339 = 0.2e1 * t1293;
t1294 = m(2) + m(3);
t1338 = g(3) * mrSges(1,2);
t1284 = cos(qJ(3,3));
t1237 = t1284 * pkin(3) + pkin(2);
t1286 = cos(qJ(3,2));
t1238 = t1286 * pkin(3) + pkin(2);
t1288 = cos(qJ(3,1));
t1239 = t1288 * pkin(3) + pkin(2);
t1275 = legFrame(3,2);
t1254 = sin(t1275);
t1257 = cos(t1275);
t1203 = t1257 * g(1) - t1254 * g(2);
t1337 = mrSges(3,2) * t1203;
t1276 = legFrame(2,2);
t1255 = sin(t1276);
t1258 = cos(t1276);
t1204 = t1258 * g(1) - t1255 * g(2);
t1336 = mrSges(3,2) * t1204;
t1277 = legFrame(1,2);
t1256 = sin(t1277);
t1259 = cos(t1277);
t1205 = t1259 * g(1) - t1256 * g(2);
t1335 = mrSges(3,2) * t1205;
t1194 = mrSges(3,1) * t1203;
t1222 = t1294 * pkin(1) + mrSges(1,1);
t1218 = g(3) * t1222;
t1221 = -0.2e1 * g(3) * t1240;
t1223 = g(3) * t1340;
t1305 = pkin(7) + qJ(3,3);
t1247 = qJ(1,3) + t1305;
t1230 = sin(t1247);
t1308 = -pkin(7) + qJ(3,3);
t1248 = qJ(1,3) - t1308;
t1231 = sin(t1248);
t1264 = qJ(1,3) + pkin(7);
t1241 = sin(t1264);
t1244 = cos(t1264);
t1290 = mrSges(3,2) * g(3);
t1291 = mrSges(3,1) * g(3);
t1185 = (t1291 - t1337) * cos(t1248) / 0.2e1 + (t1194 + t1290) * t1231 / 0.2e1 + (t1291 + t1337) * cos(t1247) / 0.2e1 + (t1194 - t1290) * t1230 / 0.2e1 + (t1203 * t1341 + t1223) * t1244 / 0.2e1 + (t1203 * t1340 + t1221) * t1241 / 0.2e1 + (mrSges(1,2) * t1203 + t1218) * t1285 + t1279 * (t1222 * t1203 - t1338);
t1334 = t1185 / 0.2e1;
t1195 = mrSges(3,1) * t1204;
t1306 = pkin(7) + qJ(3,2);
t1249 = qJ(1,2) + t1306;
t1232 = sin(t1249);
t1309 = -pkin(7) + qJ(3,2);
t1250 = qJ(1,2) - t1309;
t1233 = sin(t1250);
t1265 = qJ(1,2) + pkin(7);
t1242 = sin(t1265);
t1245 = cos(t1265);
t1186 = (t1291 - t1336) * cos(t1250) / 0.2e1 + (t1195 + t1290) * t1233 / 0.2e1 + (t1291 + t1336) * cos(t1249) / 0.2e1 + (t1195 - t1290) * t1232 / 0.2e1 + (t1204 * t1341 + t1223) * t1245 / 0.2e1 + (t1204 * t1340 + t1221) * t1242 / 0.2e1 + (mrSges(1,2) * t1204 + t1218) * t1287 + t1281 * (t1222 * t1204 - t1338);
t1333 = t1186 / 0.2e1;
t1196 = mrSges(3,1) * t1205;
t1307 = pkin(7) + qJ(3,1);
t1251 = qJ(1,1) + t1307;
t1234 = sin(t1251);
t1310 = -pkin(7) + qJ(3,1);
t1252 = qJ(1,1) - t1310;
t1235 = sin(t1252);
t1266 = qJ(1,1) + pkin(7);
t1243 = sin(t1266);
t1246 = cos(t1266);
t1187 = (t1291 - t1335) * cos(t1252) / 0.2e1 + (t1196 + t1290) * t1235 / 0.2e1 + (t1291 + t1335) * cos(t1251) / 0.2e1 + (t1196 - t1290) * t1234 / 0.2e1 + (t1205 * t1341 + t1223) * t1246 / 0.2e1 + (t1205 * t1340 + t1221) * t1243 / 0.2e1 + (mrSges(1,2) * t1205 + t1218) * t1289 + t1283 * (t1222 * t1205 - t1338);
t1332 = t1187 / 0.2e1;
t1200 = t1254 * g(1) + t1257 * g(2);
t1278 = sin(qJ(3,3));
t1295 = 0.1e1 / pkin(3);
t1331 = (-t1200 * (t1284 * mrSges(3,1) - t1278 * mrSges(3,2)) + (-g(3) * t1241 + t1203 * t1244) * (mrSges(3,1) * t1278 + mrSges(3,2) * t1284)) * t1295;
t1201 = t1255 * g(1) + t1258 * g(2);
t1280 = sin(qJ(3,2));
t1330 = (-t1201 * (t1286 * mrSges(3,1) - t1280 * mrSges(3,2)) + (-g(3) * t1242 + t1204 * t1245) * (mrSges(3,1) * t1280 + mrSges(3,2) * t1286)) * t1295;
t1202 = t1256 * g(1) + t1259 * g(2);
t1282 = sin(qJ(3,1));
t1329 = (-t1202 * (t1288 * mrSges(3,1) - t1282 * mrSges(3,2)) + (-g(3) * t1243 + t1205 * t1246) * (mrSges(3,1) * t1282 + mrSges(3,2) * t1288)) * t1295;
t1328 = t1200 * t1294;
t1327 = t1201 * t1294;
t1326 = t1202 * t1294;
t1325 = t1237 * t1285;
t1324 = t1238 * t1287;
t1323 = t1239 * t1289;
t1322 = t1278 * t1254;
t1321 = t1278 * t1257;
t1318 = t1280 * t1255;
t1317 = t1280 * t1258;
t1314 = t1282 * t1256;
t1313 = t1282 * t1259;
t1304 = pkin(3) * (-t1285 * t1274 + t1320) * t1284 ^ 2;
t1303 = pkin(3) * (-t1287 * t1274 + t1316) * t1286 ^ 2;
t1302 = pkin(3) * (-t1289 * t1274 + t1312) * t1288 ^ 2;
t1301 = ((-t1319 + t1325) * t1274 - t1206 - t1237 * t1320) * t1331;
t1300 = ((-t1315 + t1324) * t1274 - t1207 - t1238 * t1316) * t1330;
t1299 = ((-t1311 + t1323) * t1274 - t1208 - t1239 * t1312) * t1329;
t1269 = 0.1e1 / t1282;
t1268 = 0.1e1 / t1280;
t1267 = 0.1e1 / t1278;
t1253 = t1274 * pkin(1);
t1236 = t1253 + pkin(2);
t1229 = -t1277 + t1266;
t1228 = t1277 + t1266;
t1227 = -t1276 + t1265;
t1226 = t1276 + t1265;
t1225 = -t1275 + t1264;
t1224 = t1275 + t1264;
t1211 = 0.1e1 / (t1253 + t1239);
t1210 = 0.1e1 / (t1253 + t1238);
t1209 = 0.1e1 / (t1253 + t1237);
t1 = [-g(1) * m(4) + ((-sin(t1228) - sin(t1229)) * t1332 + (-(-t1259 * t1302 + (pkin(3) * t1314 - t1344 * t1259) * t1288 + t1236 * t1314) * t1326 - t1259 * t1299) * t1269) * t1211 + ((-sin(t1226) - sin(t1227)) * t1333 + (-(-t1258 * t1303 + (pkin(3) * t1318 - t1345 * t1258) * t1286 + t1236 * t1318) * t1327 - t1258 * t1300) * t1268) * t1210 + ((-sin(t1224) - sin(t1225)) * t1334 + (-(-t1257 * t1304 + (pkin(3) * t1322 - t1346 * t1257) * t1284 + t1236 * t1322) * t1328 - t1257 * t1301) * t1267) * t1209; -g(2) * m(4) + ((cos(t1229) - cos(t1228)) * t1332 + (-(t1256 * t1302 + (pkin(3) * t1313 + t1344 * t1256) * t1288 + t1236 * t1313) * t1326 + t1256 * t1299) * t1269) * t1211 + ((cos(t1227) - cos(t1226)) * t1333 + (-(t1255 * t1303 + (pkin(3) * t1317 + t1345 * t1255) * t1286 + t1236 * t1317) * t1327 + t1255 * t1300) * t1268) * t1210 + ((cos(t1225) - cos(t1224)) * t1334 + (-(t1254 * t1304 + (pkin(3) * t1321 + t1346 * t1254) * t1284 + t1236 * t1321) * t1328 + t1254 * t1301) * t1267) * t1209; (t1246 * t1339 + t1283 * t1343 + t1243 * t1342 + (t1234 + t1235) * pkin(3)) / (pkin(3) * sin(0.2e1 * qJ(3,1)) + t1282 * t1342 + (sin(t1307) + sin(t1310)) * pkin(1)) * t1329 + (t1245 * t1339 + t1281 * t1343 + t1242 * t1342 + (t1232 + t1233) * pkin(3)) / (pkin(3) * sin(0.2e1 * qJ(3,2)) + t1280 * t1342 + (sin(t1306) + sin(t1309)) * pkin(1)) * t1330 + (t1244 * t1339 + t1279 * t1343 + t1241 * t1342 + (t1230 + t1231) * pkin(3)) / (pkin(3) * sin(0.2e1 * qJ(3,3)) + t1278 * t1342 + (sin(t1305) + sin(t1308)) * pkin(1)) * t1331 - g(3) * m(4) + (-t1246 * t1187 + t1288 * ((t1239 * t1283 + t1289 * t1293) * t1274 - t1219 * t1283 + t1273 * t1323) * t1269 * t1326) * t1211 + (-t1245 * t1186 + t1286 * ((t1238 * t1281 + t1287 * t1293) * t1274 - t1219 * t1281 + t1273 * t1324) * t1268 * t1327) * t1210 + (-t1244 * t1185 + t1284 * ((t1237 * t1279 + t1285 * t1293) * t1274 - t1219 * t1279 + t1273 * t1325) * t1267 * t1328) * t1209;];
taugX  = t1;
