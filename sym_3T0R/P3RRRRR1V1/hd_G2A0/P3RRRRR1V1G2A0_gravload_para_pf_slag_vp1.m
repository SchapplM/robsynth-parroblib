% Calculate Gravitation load for parallel robot
% P3RRRRR1V1G2A0
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1]';
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
% Datum: 2020-08-07 03:36
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRRRR1V1G2A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(4,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR1V1G2A0_gravload_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR1V1G2A0_gravload_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RRRRR1V1G2A0_gravload_para_pf_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRRRR1V1G2A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRRRR1V1G2A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRRRR1V1G2A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR1V1G2A0_gravload_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR1V1G2A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 03:33:40
% EndTime: 2020-08-07 03:33:41
% DurationCPUTime: 0.84s
% Computational Cost: add. (717->134), mult. (1422->219), div. (54->8), fcn. (816->36), ass. (0->104)
t1288 = m(1) * rSges(1,2) + m(2) * rSges(2,3) + m(3) * rSges(3,3);
t1326 = m(1) * rSges(1,1) + (m(2) + m(3)) * pkin(1);
t1255 = sin(qJ(3,1));
t1264 = cos(qJ(3,1));
t1290 = -m(2) * rSges(2,1) - pkin(2) * m(3);
t1319 = m(3) * rSges(3,2);
t1320 = m(3) * rSges(3,1);
t1208 = t1255 * t1319 - t1264 * t1320 + t1290;
t1268 = m(2) * rSges(2,2);
t1211 = t1268 + (rSges(3,1) * t1255 + rSges(3,2) * t1264) * m(3);
t1256 = sin(qJ(2,1));
t1265 = cos(qJ(2,1));
t1323 = -t1208 * t1265 - t1211 * t1256;
t1252 = sin(qJ(3,2));
t1261 = cos(qJ(3,2));
t1207 = t1252 * t1319 - t1261 * t1320 + t1290;
t1210 = t1268 + (rSges(3,1) * t1252 + rSges(3,2) * t1261) * m(3);
t1253 = sin(qJ(2,2));
t1262 = cos(qJ(2,2));
t1322 = -t1207 * t1262 - t1210 * t1253;
t1249 = sin(qJ(3,3));
t1258 = cos(qJ(3,3));
t1206 = t1249 * t1319 - t1258 * t1320 + t1290;
t1209 = t1268 + (rSges(3,1) * t1249 + rSges(3,2) * t1258) * m(3);
t1250 = sin(qJ(2,3));
t1259 = cos(qJ(2,3));
t1321 = -t1206 * t1259 - t1209 * t1250;
t1318 = m(3) / pkin(3);
t1245 = qJ(2,1) + qJ(3,1);
t1244 = qJ(2,2) + qJ(3,2);
t1243 = qJ(2,3) + qJ(3,3);
t1246 = legFrame(3,2);
t1234 = sin(t1246);
t1237 = cos(t1246);
t1218 = t1237 * g(1) - t1234 * g(2);
t1231 = cos(t1243);
t1251 = sin(qJ(1,3));
t1260 = cos(qJ(1,3));
t1283 = t1326 * g(3);
t1284 = t1288 * g(3);
t1317 = (t1284 * t1260 + (t1321 * g(3) + t1283) * t1251 + ((-t1326 - t1321) * t1260 + t1288 * t1251) * t1218) / (t1259 * pkin(2) + pkin(3) * t1231 + pkin(1));
t1247 = legFrame(2,2);
t1235 = sin(t1247);
t1238 = cos(t1247);
t1219 = t1238 * g(1) - t1235 * g(2);
t1232 = cos(t1244);
t1254 = sin(qJ(1,2));
t1263 = cos(qJ(1,2));
t1316 = (t1284 * t1263 + (t1322 * g(3) + t1283) * t1254 + ((-t1326 - t1322) * t1263 + t1288 * t1254) * t1219) / (t1262 * pkin(2) + pkin(3) * t1232 + pkin(1));
t1248 = legFrame(1,2);
t1236 = sin(t1248);
t1239 = cos(t1248);
t1220 = t1239 * g(1) - t1236 * g(2);
t1233 = cos(t1245);
t1257 = sin(qJ(1,1));
t1266 = cos(qJ(1,1));
t1315 = (t1284 * t1266 + (t1323 * g(3) + t1283) * t1257 + ((-t1326 - t1323) * t1266 + t1288 * t1257) * t1220) / (t1265 * pkin(2) + pkin(3) * t1233 + pkin(1));
t1215 = t1234 * g(1) + t1237 * g(2);
t1240 = 0.1e1 / t1249;
t1282 = g(3) * t1260 + t1218 * t1251;
t1314 = ((-t1282 * t1206 - t1215 * t1209) * t1250 + t1259 * (-t1215 * t1206 + t1282 * t1209)) * t1240;
t1216 = t1235 * g(1) + t1238 * g(2);
t1241 = 0.1e1 / t1252;
t1281 = g(3) * t1263 + t1219 * t1254;
t1313 = ((-t1281 * t1207 - t1216 * t1210) * t1253 + t1262 * (-t1216 * t1207 + t1281 * t1210)) * t1241;
t1217 = t1236 * g(1) + t1239 * g(2);
t1242 = 0.1e1 / t1255;
t1280 = g(3) * t1266 + t1220 * t1257;
t1312 = ((-t1280 * t1208 - t1217 * t1211) * t1256 + t1265 * (-t1217 * t1208 + t1280 * t1211)) * t1242;
t1311 = ((rSges(3,1) * t1215 + t1282 * rSges(3,2)) * t1231 + sin(t1243) * (t1282 * rSges(3,1) - rSges(3,2) * t1215)) * t1240;
t1310 = ((rSges(3,1) * t1216 + t1281 * rSges(3,2)) * t1232 + sin(t1244) * (t1281 * rSges(3,1) - rSges(3,2) * t1216)) * t1241;
t1309 = ((rSges(3,1) * t1217 + t1280 * rSges(3,2)) * t1233 + sin(t1245) * (t1280 * rSges(3,1) - rSges(3,2) * t1217)) * t1242;
t1302 = t1234 * t1251;
t1301 = t1235 * t1254;
t1300 = t1236 * t1257;
t1299 = t1237 * t1251;
t1298 = t1238 * t1254;
t1297 = t1239 * t1257;
t1296 = t1249 * t1250;
t1295 = t1249 * t1259;
t1294 = t1252 * t1253;
t1293 = t1252 * t1262;
t1292 = t1255 * t1256;
t1291 = t1255 * t1265;
t1287 = t1260 * t1317;
t1286 = t1263 * t1316;
t1285 = t1266 * t1315;
t1279 = -t1250 * t1258 - t1295;
t1278 = -t1258 * t1259 + t1296;
t1277 = -t1253 * t1261 - t1293;
t1276 = -t1261 * t1262 + t1294;
t1275 = -t1256 * t1264 - t1291;
t1274 = -t1264 * t1265 + t1292;
t1273 = 0.1e1 / pkin(2);
t1227 = pkin(3) * t1264 + pkin(2);
t1226 = pkin(3) * t1261 + pkin(2);
t1225 = pkin(3) * t1258 + pkin(2);
t1205 = -pkin(3) * t1292 + t1227 * t1265;
t1204 = -pkin(3) * t1294 + t1226 * t1262;
t1203 = -pkin(3) * t1296 + t1225 * t1259;
t1202 = pkin(3) * t1291 + t1256 * t1227;
t1201 = pkin(3) * t1293 + t1253 * t1226;
t1200 = pkin(3) * t1295 + t1250 * t1225;
t1 = [t1237 * t1287 + t1238 * t1286 + t1239 * t1285 - m(4) * g(1) + ((t1275 * t1236 - t1274 * t1297) * t1312 + (t1277 * t1235 - t1276 * t1298) * t1313 + (t1279 * t1234 - t1278 * t1299) * t1314 + ((t1236 * t1202 - t1205 * t1297) * t1309 + (t1235 * t1201 - t1204 * t1298) * t1310 + (t1234 * t1200 - t1203 * t1299) * t1311) * t1318) * t1273; -t1234 * t1287 - t1235 * t1286 - t1236 * t1285 - m(4) * g(2) + ((t1275 * t1239 + t1274 * t1300) * t1312 + (t1277 * t1238 + t1276 * t1301) * t1313 + (t1279 * t1237 + t1278 * t1302) * t1314 + ((t1202 * t1239 + t1205 * t1300) * t1309 + (t1201 * t1238 + t1204 * t1301) * t1310 + (t1200 * t1237 + t1203 * t1302) * t1311) * t1318) * t1273; -t1251 * t1317 - t1254 * t1316 - t1257 * t1315 - m(4) * g(3) + ((cos(qJ(1,1) - t1245) + cos(qJ(1,1) + t1245)) * t1312 / 0.2e1 + (cos(qJ(1,2) - t1244) + cos(qJ(1,2) + t1244)) * t1313 / 0.2e1 + (cos(qJ(1,3) - t1243) + cos(qJ(1,3) + t1243)) * t1314 / 0.2e1 + (-t1203 * t1260 * t1311 - t1204 * t1263 * t1310 - t1205 * t1266 * t1309) * t1318) * t1273;];
taugX  = t1;
