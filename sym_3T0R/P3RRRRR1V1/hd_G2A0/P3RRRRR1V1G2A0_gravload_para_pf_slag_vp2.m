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
% Datum: 2020-08-07 03:36
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRRRR1V1G2A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(4,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR1V1G2A0_gravload_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR1V1G2A0_gravload_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RRRRR1V1G2A0_gravload_para_pf_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRRRR1V1G2A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRRRR1V1G2A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRRRR1V1G2A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR1V1G2A0_gravload_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR1V1G2A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 03:34:06
% EndTime: 2020-08-07 03:34:06
% DurationCPUTime: 0.63s
% Computational Cost: add. (699->133), mult. (1017->208), div. (54->8), fcn. (780->36), ass. (0->101)
t1229 = sin(qJ(3,1));
t1238 = cos(qJ(3,1));
t1283 = -m(3) * pkin(2) - mrSges(2,1);
t1182 = -mrSges(3,1) * t1238 + mrSges(3,2) * t1229 + t1283;
t1194 = t1229 * mrSges(3,1) + t1238 * mrSges(3,2) + mrSges(2,2);
t1230 = sin(qJ(2,1));
t1239 = cos(qJ(2,1));
t1286 = -t1239 * t1182 - t1194 * t1230;
t1226 = sin(qJ(3,2));
t1235 = cos(qJ(3,2));
t1181 = -mrSges(3,1) * t1235 + mrSges(3,2) * t1226 + t1283;
t1193 = t1226 * mrSges(3,1) + t1235 * mrSges(3,2) + mrSges(2,2);
t1227 = sin(qJ(2,2));
t1236 = cos(qJ(2,2));
t1285 = -t1236 * t1181 - t1193 * t1227;
t1223 = sin(qJ(3,3));
t1232 = cos(qJ(3,3));
t1180 = -mrSges(3,1) * t1232 + mrSges(3,2) * t1223 + t1283;
t1192 = t1223 * mrSges(3,1) + t1232 * mrSges(3,2) + mrSges(2,2);
t1224 = sin(qJ(2,3));
t1233 = cos(qJ(2,3));
t1284 = -t1233 * t1180 - t1192 * t1224;
t1219 = qJ(2,1) + qJ(3,1);
t1218 = qJ(2,2) + qJ(3,2);
t1217 = qJ(2,3) + qJ(3,3);
t1220 = legFrame(3,2);
t1207 = sin(t1220);
t1210 = cos(t1220);
t1189 = t1210 * g(1) - t1207 * g(2);
t1196 = mrSges(1,1) + (m(2) + m(3)) * pkin(1);
t1195 = t1196 * g(3);
t1213 = mrSges(2,3) + mrSges(3,3) + mrSges(1,2);
t1200 = t1213 * g(3);
t1201 = cos(t1217);
t1225 = sin(qJ(1,3));
t1234 = cos(qJ(1,3));
t1282 = (t1200 * t1234 + (t1284 * g(3) + t1195) * t1225 + ((-t1196 - t1284) * t1234 + t1213 * t1225) * t1189) / (t1233 * pkin(2) + pkin(3) * t1201 + pkin(1));
t1221 = legFrame(2,2);
t1208 = sin(t1221);
t1211 = cos(t1221);
t1190 = t1211 * g(1) - t1208 * g(2);
t1202 = cos(t1218);
t1228 = sin(qJ(1,2));
t1237 = cos(qJ(1,2));
t1281 = (t1200 * t1237 + (t1285 * g(3) + t1195) * t1228 + ((-t1196 - t1285) * t1237 + t1213 * t1228) * t1190) / (t1236 * pkin(2) + pkin(3) * t1202 + pkin(1));
t1222 = legFrame(1,2);
t1209 = sin(t1222);
t1212 = cos(t1222);
t1191 = t1212 * g(1) - t1209 * g(2);
t1203 = cos(t1219);
t1231 = sin(qJ(1,1));
t1240 = cos(qJ(1,1));
t1280 = (t1200 * t1240 + (t1286 * g(3) + t1195) * t1231 + ((-t1196 - t1286) * t1240 + t1213 * t1231) * t1191) / (t1239 * pkin(2) + pkin(3) * t1203 + pkin(1));
t1186 = t1207 * g(1) + t1210 * g(2);
t1214 = 0.1e1 / t1223;
t1252 = g(3) * t1234 + t1189 * t1225;
t1279 = ((-t1252 * t1180 - t1186 * t1192) * t1224 + t1233 * (-t1186 * t1180 + t1252 * t1192)) * t1214;
t1187 = t1208 * g(1) + t1211 * g(2);
t1215 = 0.1e1 / t1226;
t1251 = g(3) * t1237 + t1190 * t1228;
t1278 = ((-t1251 * t1181 - t1187 * t1193) * t1227 + t1236 * (-t1187 * t1181 + t1251 * t1193)) * t1215;
t1188 = t1209 * g(1) + t1212 * g(2);
t1216 = 0.1e1 / t1229;
t1250 = g(3) * t1240 + t1191 * t1231;
t1277 = ((-t1250 * t1182 - t1188 * t1194) * t1230 + t1239 * (-t1188 * t1182 + t1250 * t1194)) * t1216;
t1276 = ((mrSges(3,1) * t1186 + t1252 * mrSges(3,2)) * t1201 + sin(t1217) * (t1252 * mrSges(3,1) - mrSges(3,2) * t1186)) * t1214;
t1275 = ((mrSges(3,1) * t1187 + t1251 * mrSges(3,2)) * t1202 + sin(t1218) * (t1251 * mrSges(3,1) - mrSges(3,2) * t1187)) * t1215;
t1274 = ((mrSges(3,1) * t1188 + t1250 * mrSges(3,2)) * t1203 + sin(t1219) * (t1250 * mrSges(3,1) - mrSges(3,2) * t1188)) * t1216;
t1270 = t1207 * t1225;
t1269 = t1208 * t1228;
t1268 = t1209 * t1231;
t1267 = t1210 * t1225;
t1266 = t1211 * t1228;
t1265 = t1212 * t1231;
t1264 = t1223 * t1224;
t1263 = t1223 * t1233;
t1262 = t1226 * t1227;
t1261 = t1226 * t1236;
t1260 = t1229 * t1230;
t1259 = t1229 * t1239;
t1255 = t1234 * t1282;
t1254 = t1237 * t1281;
t1253 = t1240 * t1280;
t1249 = -t1224 * t1232 - t1263;
t1248 = -t1232 * t1233 + t1264;
t1247 = -t1227 * t1235 - t1261;
t1246 = -t1235 * t1236 + t1262;
t1245 = -t1230 * t1238 - t1259;
t1244 = -t1238 * t1239 + t1260;
t1243 = 0.1e1 / pkin(2);
t1242 = 0.1e1 / pkin(3);
t1199 = pkin(3) * t1238 + pkin(2);
t1198 = pkin(3) * t1235 + pkin(2);
t1197 = pkin(3) * t1232 + pkin(2);
t1179 = -pkin(3) * t1260 + t1199 * t1239;
t1178 = -pkin(3) * t1262 + t1198 * t1236;
t1177 = -pkin(3) * t1264 + t1197 * t1233;
t1176 = pkin(3) * t1259 + t1230 * t1199;
t1175 = pkin(3) * t1261 + t1227 * t1198;
t1174 = pkin(3) * t1263 + t1224 * t1197;
t1 = [t1210 * t1255 + t1211 * t1254 + t1212 * t1253 - g(1) * m(4) + ((t1245 * t1209 - t1244 * t1265) * t1277 + (t1247 * t1208 - t1246 * t1266) * t1278 + (t1249 * t1207 - t1248 * t1267) * t1279 + ((t1209 * t1176 - t1179 * t1265) * t1274 + (t1208 * t1175 - t1178 * t1266) * t1275 + (t1207 * t1174 - t1177 * t1267) * t1276) * t1242) * t1243; -t1207 * t1255 - t1208 * t1254 - t1209 * t1253 - g(2) * m(4) + ((t1245 * t1212 + t1244 * t1268) * t1277 + (t1247 * t1211 + t1246 * t1269) * t1278 + (t1249 * t1210 + t1248 * t1270) * t1279 + ((t1176 * t1212 + t1179 * t1268) * t1274 + (t1175 * t1211 + t1178 * t1269) * t1275 + (t1174 * t1210 + t1177 * t1270) * t1276) * t1242) * t1243; -t1225 * t1282 - t1228 * t1281 - t1231 * t1280 - g(3) * m(4) + ((cos(qJ(1,1) - t1219) + cos(qJ(1,1) + t1219)) * t1277 / 0.2e1 + (cos(qJ(1,2) - t1218) + cos(qJ(1,2) + t1218)) * t1278 / 0.2e1 + (cos(qJ(1,3) - t1217) + cos(qJ(1,3) + t1217)) * t1279 / 0.2e1 + (-t1234 * t1177 * t1276 - t1237 * t1178 * t1275 - t1240 * t1179 * t1274) * t1242) * t1243;];
taugX  = t1;
