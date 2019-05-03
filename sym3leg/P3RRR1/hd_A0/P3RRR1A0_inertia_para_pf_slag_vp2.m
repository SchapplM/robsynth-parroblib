% Calculate inertia matrix for parallel robot
% P3RRR1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% qJ [2x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2,d3]';
% m [3x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [3x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Ifges [3x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% MX [3x3]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 15:38
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RRR1A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(2,3),zeros(3,3),zeros(3,3),zeros(5,1),zeros(2+1,1),zeros(2+1,3),zeros(2+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRR1A0_inertia_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3RRR1A0_inertia_para_pf_slag_vp2: qJ has to be [2x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRR1A0_inertia_para_pf_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'P3RRR1A0_inertia_para_pf_slag_vp2: m has to be [3x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [3,3]), ...
  'P3RRR1A0_inertia_para_pf_slag_vp2: mrSges has to be [3x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [3 6]), ...
  'P3RRR1A0_inertia_para_pf_slag_vp2: Ifges has to be [3x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRR1A0_inertia_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRR1A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 15:38:20
% EndTime: 2019-05-03 15:38:20
% DurationCPUTime: 0.47s
% Computational Cost: add. (2000->119), mult. (2330->238), div. (540->5), fcn. (2546->26), ass. (0->123)
t310 = qJ(1,3) + qJ(2,3);
t296 = sin(t310);
t299 = cos(t310);
t317 = sin(qJ(1,3));
t323 = cos(qJ(1,3));
t376 = 0.1e1 / (t296 * t323 - t317 * t299);
t311 = qJ(1,2) + qJ(2,2);
t297 = sin(t311);
t300 = cos(t311);
t319 = sin(qJ(1,2));
t325 = cos(qJ(1,2));
t375 = 0.1e1 / (t297 * t325 - t319 * t300);
t312 = qJ(1,1) + qJ(2,1);
t298 = sin(t312);
t301 = cos(t312);
t321 = sin(qJ(1,1));
t327 = cos(qJ(1,1));
t374 = 0.1e1 / (t298 * t327 - t321 * t301);
t313 = legFrame(3,3);
t302 = sin(t313);
t305 = cos(t313);
t268 = t305 * t296 + t302 * t299;
t262 = pkin(1) * (t302 * t323 + t305 * t317) + t268 * pkin(2);
t373 = t262 * t376;
t314 = legFrame(2,3);
t303 = sin(t314);
t306 = cos(t314);
t270 = t306 * t297 + t303 * t300;
t263 = pkin(1) * (t303 * t325 + t306 * t319) + t270 * pkin(2);
t372 = t263 * t375;
t315 = legFrame(1,3);
t304 = sin(t315);
t307 = cos(t315);
t272 = t307 * t298 + t304 * t301;
t264 = pkin(1) * (t304 * t327 + t307 * t321) + t272 * pkin(2);
t371 = t264 * t374;
t269 = -t296 * t302 + t305 * t299;
t265 = pkin(1) * (-t317 * t302 + t305 * t323) + t269 * pkin(2);
t370 = t265 * t376;
t271 = -t297 * t303 + t306 * t300;
t266 = pkin(1) * (-t319 * t303 + t306 * t325) + t271 * pkin(2);
t369 = t266 * t375;
t273 = -t298 * t304 + t307 * t301;
t267 = pkin(1) * (-t321 * t304 + t307 * t327) + t273 * pkin(2);
t368 = t267 * t374;
t367 = t268 * t376;
t366 = t269 * t376;
t365 = t270 * t375;
t364 = t271 * t375;
t363 = t272 * t374;
t362 = t273 * t374;
t342 = (mrSges(2,1) * cos(qJ(2,3)) - mrSges(2,2) * sin(qJ(2,3))) * pkin(1);
t351 = m(2) * pkin(1) ^ 2 + Ifges(1,3) + Ifges(2,3);
t283 = 0.2e1 * t342 + t351;
t361 = t376 * t283;
t292 = Ifges(2,3) + t342;
t360 = t376 * t292;
t341 = (mrSges(2,1) * cos(qJ(2,2)) - mrSges(2,2) * sin(qJ(2,2))) * pkin(1);
t284 = 0.2e1 * t341 + t351;
t359 = t375 * t284;
t293 = Ifges(2,3) + t341;
t358 = t375 * t293;
t340 = (mrSges(2,1) * cos(qJ(2,1)) - mrSges(2,2) * sin(qJ(2,1))) * pkin(1);
t285 = 0.2e1 * t340 + t351;
t357 = t374 * t285;
t294 = Ifges(2,3) + t340;
t356 = t374 * t294;
t338 = 0.1e1 / pkin(2);
t355 = t376 * t338;
t354 = t375 * t338;
t353 = t374 * t338;
t339 = 0.1e1 / pkin(1);
t352 = t338 * t339;
t350 = Ifges(2,3) * t355;
t349 = Ifges(2,3) * t354;
t348 = Ifges(2,3) * t353;
t347 = t292 * t355;
t346 = t293 * t354;
t345 = t294 * t353;
t329 = xP(3);
t308 = sin(t329);
t309 = cos(t329);
t330 = mrSges(3,2);
t331 = mrSges(3,1);
t344 = -t308 * t330 + t309 * t331;
t343 = -t308 * t331 - t309 * t330;
t337 = koppelP(1,1);
t336 = koppelP(2,1);
t335 = koppelP(3,1);
t334 = koppelP(1,2);
t333 = koppelP(2,2);
t332 = koppelP(3,2);
t291 = -t308 * t334 + t309 * t337;
t290 = -t308 * t333 + t309 * t336;
t289 = -t308 * t332 + t309 * t335;
t288 = -t308 * t337 - t309 * t334;
t287 = -t308 * t336 - t309 * t333;
t286 = -t308 * t335 - t309 * t332;
t261 = (t272 * t291 + t273 * t288) * t339 * t374;
t260 = (t270 * t290 + t271 * t287) * t339 * t375;
t259 = (t268 * t289 + t269 * t286) * t339 * t376;
t258 = (-t267 * t348 + t273 * t356) * t339;
t257 = (-t264 * t348 + t272 * t356) * t339;
t256 = (-t266 * t349 + t271 * t358) * t339;
t255 = (-t263 * t349 + t270 * t358) * t339;
t254 = (-t265 * t350 + t269 * t360) * t339;
t253 = (-t262 * t350 + t268 * t360) * t339;
t252 = (-t267 * t345 + t273 * t357) * t339;
t251 = (-t264 * t345 + t272 * t357) * t339;
t250 = (-t266 * t346 + t271 * t359) * t339;
t249 = (-t263 * t346 + t270 * t359) * t339;
t248 = (-t265 * t347 + t269 * t361) * t339;
t247 = (-t262 * t347 + t268 * t361) * t339;
t246 = (t264 * t291 + t267 * t288) * t374 * t352;
t245 = (t263 * t290 + t266 * t287) * t375 * t352;
t244 = (t262 * t289 + t265 * t286) * t376 * t352;
t243 = -t246 * Ifges(2,3) + t261 * t294;
t242 = -t245 * Ifges(2,3) + t260 * t293;
t241 = -t244 * Ifges(2,3) + t259 * t292;
t240 = -t246 * t294 + t261 * t285;
t239 = -t245 * t293 + t260 * t284;
t238 = -t244 * t292 + t259 * t283;
t1 = [m(3) + (t248 * t366 + t250 * t364 + t252 * t362 + (-t254 * t370 - t256 * t369 - t258 * t368) * t338) * t339, (t248 * t367 + t250 * t365 + t252 * t363 + (-t254 * t373 - t256 * t372 - t258 * t371) * t338) * t339, -t254 * t244 - t256 * t245 - t258 * t246 + t248 * t259 + t250 * t260 + t252 * t261 + t343; (t247 * t366 + t249 * t364 + t251 * t362 + (-t253 * t370 - t255 * t369 - t257 * t368) * t338) * t339, m(3) + (t247 * t367 + t249 * t365 + t251 * t363 + (-t253 * t373 - t255 * t372 - t257 * t371) * t338) * t339, -t253 * t244 - t255 * t245 - t257 * t246 + t247 * t259 + t249 * t260 + t251 * t261 + t344; (t238 * t366 + t239 * t364 + t240 * t362 + (-t241 * t370 - t242 * t369 - t243 * t368) * t338) * t339 + t343, (t238 * t367 + t239 * t365 + t240 * t363 + (-t241 * t373 - t242 * t372 - t243 * t371) * t338) * t339 + t344, t238 * t259 + t239 * t260 + t240 * t261 - t241 * t244 - t242 * t245 - t243 * t246 + Ifges(3,3);];
MX  = t1;
