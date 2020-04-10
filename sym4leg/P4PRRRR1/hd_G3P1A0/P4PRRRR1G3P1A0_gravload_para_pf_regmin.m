% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P4PRRRR1G3P1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [4x1]
%   Generalized platform coordinates
% qJ [3x4]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [4x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% koppelP [4x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [4x15]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-02 19:06
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P4PRRRR1G3P1A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(3,1),zeros(4,3),zeros(4,3),zeros(2,1)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR1G3P1A0_gravload_para_pf_regmin: xP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR1G3P1A0_gravload_para_pf_regmin: qJ has to be [3x4] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P4PRRRR1G3P1A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P4PRRRR1G3P1A0_gravload_para_pf_regmin: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR1G3P1A0_gravload_para_pf_regmin: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR1G3P1A0_gravload_para_pf_regmin: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-02 19:05:33
% EndTime: 2020-03-02 19:05:34
% DurationCPUTime: 1.46s
% Computational Cost: add. (393->133), mult. (804->284), div. (220->13), fcn. (932->26), ass. (0->135)
t1248 = legFrame(4,2);
t1222 = sin(t1248);
t1226 = cos(t1248);
t1212 = t1222 * g(1) + t1226 * g(2);
t1216 = t1226 * g(1) - t1222 * g(2);
t1245 = sin(qJ(2,4));
t1247 = cos(qJ(2,4));
t1190 = t1212 * t1247 - t1216 * t1245;
t1232 = 0.1e1 / t1245;
t1325 = t1190 * t1232;
t1249 = legFrame(3,2);
t1223 = sin(t1249);
t1227 = cos(t1249);
t1213 = t1223 * g(1) + t1227 * g(2);
t1217 = t1227 * g(1) - t1223 * g(2);
t1253 = sin(qJ(2,3));
t1259 = cos(qJ(2,3));
t1193 = t1213 * t1259 - t1217 * t1253;
t1235 = 0.1e1 / t1253;
t1324 = t1193 * t1235;
t1250 = legFrame(2,2);
t1224 = sin(t1250);
t1228 = cos(t1250);
t1214 = t1224 * g(1) + t1228 * g(2);
t1218 = t1228 * g(1) - t1224 * g(2);
t1255 = sin(qJ(2,2));
t1261 = cos(qJ(2,2));
t1196 = t1214 * t1261 - t1218 * t1255;
t1236 = 0.1e1 / t1255;
t1323 = t1196 * t1236;
t1251 = legFrame(1,2);
t1225 = sin(t1251);
t1229 = cos(t1251);
t1215 = t1225 * g(1) + t1229 * g(2);
t1219 = t1229 * g(1) - t1225 * g(2);
t1257 = sin(qJ(2,1));
t1263 = cos(qJ(2,1));
t1199 = t1215 * t1263 - t1219 * t1257;
t1237 = 0.1e1 / t1257;
t1322 = t1199 * t1237;
t1321 = t1212 * t1232;
t1320 = t1213 * t1235;
t1319 = t1214 * t1236;
t1318 = t1215 * t1237;
t1246 = cos(qJ(3,4));
t1233 = 0.1e1 / t1246;
t1317 = t1232 * t1233;
t1316 = t1232 * t1247;
t1258 = cos(qJ(3,3));
t1238 = 0.1e1 / t1258;
t1315 = t1235 * t1238;
t1314 = t1235 * t1259;
t1260 = cos(qJ(3,2));
t1240 = 0.1e1 / t1260;
t1313 = t1236 * t1240;
t1312 = t1236 * t1261;
t1262 = cos(qJ(3,1));
t1242 = 0.1e1 / t1262;
t1311 = t1237 * t1242;
t1310 = t1237 * t1263;
t1264 = xP(4);
t1230 = sin(t1264);
t1231 = cos(t1264);
t1265 = koppelP(4,2);
t1269 = koppelP(4,1);
t1280 = -t1230 * t1265 + t1231 * t1269;
t1281 = t1230 * t1269 + t1231 * t1265;
t1184 = t1280 * t1222 + t1281 * t1226;
t1309 = t1184 * t1317;
t1266 = koppelP(3,2);
t1270 = koppelP(3,1);
t1278 = -t1230 * t1266 + t1231 * t1270;
t1279 = t1230 * t1270 + t1231 * t1266;
t1185 = t1278 * t1223 + t1279 * t1227;
t1308 = t1185 * t1315;
t1267 = koppelP(2,2);
t1271 = koppelP(2,1);
t1276 = -t1230 * t1267 + t1231 * t1271;
t1277 = t1230 * t1271 + t1231 * t1267;
t1186 = t1276 * t1224 + t1277 * t1228;
t1307 = t1186 * t1313;
t1268 = koppelP(1,2);
t1272 = koppelP(1,1);
t1274 = -t1230 * t1268 + t1231 * t1272;
t1275 = t1230 * t1272 + t1231 * t1268;
t1187 = t1274 * t1225 + t1275 * t1229;
t1306 = t1187 * t1311;
t1305 = t1222 * t1317;
t1304 = t1223 * t1315;
t1303 = t1224 * t1313;
t1302 = t1225 * t1311;
t1301 = t1226 * t1317;
t1300 = t1227 * t1315;
t1299 = t1228 * t1313;
t1298 = t1229 * t1311;
t1244 = sin(qJ(3,4));
t1297 = t1244 * t1317;
t1296 = 0.1e1 / t1246 ^ 2 * t1316;
t1252 = sin(qJ(3,3));
t1295 = t1252 * t1315;
t1294 = 0.1e1 / t1258 ^ 2 * t1314;
t1254 = sin(qJ(3,2));
t1293 = t1254 * t1313;
t1292 = 0.1e1 / t1260 ^ 2 * t1312;
t1256 = sin(qJ(3,1));
t1291 = t1256 * t1311;
t1290 = 0.1e1 / t1262 ^ 2 * t1310;
t1289 = t1190 * t1297;
t1288 = t1193 * t1295;
t1287 = t1196 * t1293;
t1286 = t1199 * t1291;
t1285 = t1244 * t1296;
t1284 = t1252 * t1294;
t1283 = t1254 * t1292;
t1282 = t1256 * t1290;
t1273 = 0.1e1 / pkin(2);
t1221 = t1231 * g(1) + t1230 * g(2);
t1220 = t1230 * g(1) - t1231 * g(2);
t1211 = t1257 * t1225 + t1229 * t1263;
t1210 = -t1225 * t1263 + t1229 * t1257;
t1209 = t1255 * t1224 + t1228 * t1261;
t1208 = -t1224 * t1261 + t1228 * t1255;
t1207 = t1253 * t1223 + t1227 * t1259;
t1206 = -t1223 * t1259 + t1227 * t1253;
t1205 = t1245 * t1222 + t1226 * t1247;
t1204 = -t1222 * t1247 + t1226 * t1245;
t1203 = (g(1) * t1263 + g(2) * t1257) * t1229 + t1225 * (g(1) * t1257 - g(2) * t1263);
t1202 = (g(1) * t1261 + t1255 * g(2)) * t1228 + t1224 * (g(1) * t1255 - g(2) * t1261);
t1201 = (g(1) * t1259 + g(2) * t1253) * t1227 + t1223 * (g(1) * t1253 - g(2) * t1259);
t1200 = (t1247 * g(1) + t1245 * g(2)) * t1226 + t1222 * (g(1) * t1245 - g(2) * t1247);
t1198 = t1215 * t1257 + t1219 * t1263;
t1195 = t1214 * t1255 + t1218 * t1261;
t1192 = t1213 * t1253 + t1217 * t1259;
t1189 = t1212 * t1245 + t1216 * t1247;
t1 = [-t1205 * t1321 - t1207 * t1320 - t1209 * t1319 - t1211 * t1318, 0, (t1190 * t1301 + t1193 * t1300 + t1196 * t1299 + t1199 * t1298) * t1273, (-t1189 * t1301 - t1192 * t1300 - t1195 * t1299 - t1198 * t1298) * t1273, 0, 0, 0, 0, 0, (t1226 * t1325 + t1227 * t1324 + t1228 * t1323 + t1229 * t1322) * t1273, (-t1226 * t1289 - t1227 * t1288 - t1228 * t1287 - t1229 * t1286) * t1273, 0, 0, 0, -t1230 * t1220 - t1231 * t1221; -t1204 * t1321 - t1206 * t1320 - t1208 * t1319 - t1210 * t1318, 0, (-t1190 * t1305 - t1193 * t1304 - t1196 * t1303 - t1199 * t1302) * t1273, (t1189 * t1305 + t1192 * t1304 + t1195 * t1303 + t1198 * t1302) * t1273, 0, 0, 0, 0, 0, (-t1222 * t1325 - t1223 * t1324 - t1224 * t1323 - t1225 * t1322) * t1273, (t1222 * t1289 + t1223 * t1288 + t1224 * t1287 + t1225 * t1286) * t1273, 0, 0, 0, t1231 * t1220 - t1230 * t1221; -t1212 * t1297 - t1213 * t1295 - t1214 * t1293 - t1215 * t1291, 0, (t1190 * t1285 + t1193 * t1284 + t1196 * t1283 + t1199 * t1282) * t1273, (-t1189 * t1285 - t1192 * t1284 - t1195 * t1283 - t1198 * t1282) * t1273, 0, 0, 0, 0, 0, ((-g(3) * t1262 + (t1199 * t1310 + t1203) * t1256) * t1242 + (-g(3) * t1260 + (t1196 * t1312 + t1202) * t1254) * t1240 + (-g(3) * t1258 + (t1193 * t1314 + t1201) * t1252) * t1238 + (-g(3) * t1246 + (t1190 * t1316 + t1200) * t1244) * t1233) * t1273, (-t1256 ^ 2 * t1199 * t1290 + t1242 * (g(3) * t1256 + t1203 * t1262) - t1254 ^ 2 * t1196 * t1292 + t1240 * (g(3) * t1254 + t1202 * t1260) - t1252 ^ 2 * t1193 * t1294 + t1238 * (g(3) * t1252 + t1201 * t1258) - t1244 ^ 2 * t1190 * t1296 + t1233 * (g(3) * t1244 + t1200 * t1246)) * t1273, 0, 0, 0, -g(3); -(t1210 * t1274 - t1211 * t1275) * t1318 - (t1208 * t1276 - t1209 * t1277) * t1319 - (t1206 * t1278 - t1207 * t1279) * t1320 - (t1204 * t1280 - t1205 * t1281) * t1321, 0, (-t1190 * t1309 - t1193 * t1308 - t1196 * t1307 - t1199 * t1306) * t1273, (t1189 * t1309 + t1192 * t1308 + t1195 * t1307 + t1198 * t1306) * t1273, 0, 0, 0, 0, 0, (-t1184 * t1325 - t1185 * t1324 - t1186 * t1323 - t1187 * t1322) * t1273, (t1184 * t1289 + t1185 * t1288 + t1186 * t1287 + t1187 * t1286) * t1273, 0, t1220, t1221, 0;];
tau_reg  = t1;
