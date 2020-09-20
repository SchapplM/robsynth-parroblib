% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3PRRR1G1P1A0
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
%   see P3PRRR1G1P1A0_convert_par2_MPV_fixb.m

% Output:
% taugX [3x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:15
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3PRRR1G1P1A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(8,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR1G1P1A0_gravload_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR1G1P1A0_gravload_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRR1G1P1A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3PRRR1G1P1A0_gravload_para_pf_mdp: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR1G1P1A0_gravload_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR1G1P1A0_gravload_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [8 1]), ...
  'P3PRRR1G1P1A0_gravload_para_pf_mdp: MDP has to be [8x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:15:00
% EndTime: 2020-03-09 21:15:00
% DurationCPUTime: 0.29s
% Computational Cost: add. (745->82), mult. (469->146), div. (84->5), fcn. (552->18), ass. (0->75)
t1273 = legFrame(3,3);
t1264 = sin(t1273);
t1267 = cos(t1273);
t1240 = t1264 * g(1) - t1267 * g(2);
t1243 = t1267 * g(1) + t1264 * g(2);
t1270 = pkin(7) + qJ(2,3);
t1261 = qJ(3,3) + t1270;
t1246 = sin(t1261);
t1249 = cos(t1261);
t1219 = t1240 * t1249 + t1243 * t1246;
t1255 = sin(t1270);
t1258 = cos(t1270);
t1231 = 0.1e1 / (-t1246 * t1258 + t1255 * t1249);
t1292 = t1219 * t1231;
t1274 = legFrame(2,3);
t1265 = sin(t1274);
t1268 = cos(t1274);
t1241 = t1265 * g(1) - t1268 * g(2);
t1244 = t1268 * g(1) + t1265 * g(2);
t1271 = pkin(7) + qJ(2,2);
t1262 = qJ(3,2) + t1271;
t1247 = sin(t1262);
t1250 = cos(t1262);
t1220 = t1241 * t1250 + t1244 * t1247;
t1256 = sin(t1271);
t1259 = cos(t1271);
t1232 = 0.1e1 / (-t1247 * t1259 + t1256 * t1250);
t1291 = t1220 * t1232;
t1275 = legFrame(1,3);
t1266 = sin(t1275);
t1269 = cos(t1275);
t1242 = t1266 * g(1) - t1269 * g(2);
t1245 = t1269 * g(1) + t1266 * g(2);
t1272 = pkin(7) + qJ(2,1);
t1263 = qJ(3,1) + t1272;
t1248 = sin(t1263);
t1251 = cos(t1263);
t1221 = t1242 * t1251 + t1245 * t1248;
t1257 = sin(t1272);
t1260 = cos(t1272);
t1233 = 0.1e1 / (-t1248 * t1260 + t1257 * t1251);
t1290 = t1221 * t1233;
t1222 = -t1240 * t1246 + t1243 * t1249;
t1289 = t1222 * t1231;
t1223 = -t1241 * t1247 + t1244 * t1250;
t1288 = t1223 * t1232;
t1224 = -t1242 * t1248 + t1245 * t1251;
t1287 = t1224 * t1233;
t1280 = t1267 * t1246 + t1264 * t1249;
t1286 = t1231 * t1280;
t1235 = t1246 * t1264 - t1267 * t1249;
t1285 = t1231 * t1235;
t1279 = t1268 * t1247 + t1265 * t1250;
t1284 = t1232 * t1279;
t1237 = t1247 * t1265 - t1268 * t1250;
t1283 = t1232 * t1237;
t1278 = t1269 * t1248 + t1266 * t1251;
t1282 = t1233 * t1278;
t1239 = t1248 * t1266 - t1269 * t1251;
t1281 = t1233 * t1239;
t1277 = 0.1e1 / pkin(2);
t1276 = 0.1e1 / pkin(3);
t1230 = -t1242 * t1257 + t1245 * t1260;
t1229 = -t1241 * t1256 + t1244 * t1259;
t1228 = -t1240 * t1255 + t1243 * t1258;
t1227 = t1242 * t1260 + t1245 * t1257;
t1226 = t1241 * t1259 + t1244 * t1256;
t1225 = t1240 * t1258 + t1243 * t1255;
t1218 = -pkin(2) * (t1257 * t1266 - t1269 * t1260) - t1239 * pkin(3);
t1217 = -pkin(2) * (t1256 * t1265 - t1268 * t1259) - t1237 * pkin(3);
t1216 = -pkin(2) * (t1255 * t1264 - t1267 * t1258) - t1235 * pkin(3);
t1215 = pkin(2) * (t1257 * t1269 + t1266 * t1260) + t1278 * pkin(3);
t1214 = pkin(2) * (t1256 * t1268 + t1265 * t1259) + t1279 * pkin(3);
t1213 = pkin(2) * (t1255 * t1267 + t1264 * t1258) + t1280 * pkin(3);
t1 = [-g(1) * MDP(8) + ((t1225 * t1285 + t1226 * t1283 + t1227 * t1281) * MDP(3) + (t1228 * t1285 + t1229 * t1283 + t1230 * t1281) * MDP(4) + (t1219 * t1285 + t1220 * t1283 + t1221 * t1281) * MDP(6) + (t1222 * t1285 + t1223 * t1283 + t1224 * t1281) * MDP(7) + ((t1216 * t1292 + t1217 * t1291 + t1218 * t1290) * MDP(6) + (t1216 * t1289 + t1217 * t1288 + t1218 * t1287) * MDP(7)) * t1276) * t1277; -g(2) * MDP(8) + ((-t1225 * t1286 - t1226 * t1284 - t1227 * t1282) * MDP(3) + (-t1228 * t1286 - t1229 * t1284 - t1230 * t1282) * MDP(4) + (-t1219 * t1286 - t1220 * t1284 - t1221 * t1282) * MDP(6) + (-t1222 * t1286 - t1223 * t1284 - t1224 * t1282) * MDP(7) + ((t1213 * t1292 + t1214 * t1291 + t1215 * t1290) * MDP(6) + (t1213 * t1289 + t1214 * t1288 + t1215 * t1287) * MDP(7)) * t1276) * t1277; (-(3 * MDP(1)) - MDP(8)) * g(3);];
taugX  = t1;
