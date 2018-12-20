% Calculate vector of centrifugal and coriolis load on the joints for
% P3RPR1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d3]';
% m [3x1]
%   mass of all robot links (including platform)
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
% taucX [3x1]
%   forces required to compensate Coriolis and centrifugal joint torques
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2018-12-20 17:54
% Revision: f9720dcdc4676342702b46a014e894344751412a
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taucX = P3RPR1A0_coriolisvec_para_pf_slag_vp2(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(2,3),zeros(3,3),zeros(3,3),zeros(4,1),zeros(2+1,1),zeros(2+1,3),zeros(2+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPR1A0_coriolisvec_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPR1A0_coriolisvec_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3RPR1A0_coriolisvec_para_pf_slag_vp2: qJ has to be [2x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RPR1A0_coriolisvec_para_pf_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'P3RPR1A0_coriolisvec_para_pf_slag_vp2: m has to be [3x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [3,3]), ...
  'P3RPR1A0_coriolisvec_para_pf_slag_vp2: mrSges has to be [3x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [3 6]), ...
  'P3RPR1A0_coriolisvec_para_pf_slag_vp2: Ifges has to be [3x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPR1A0_coriolisvec_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPR1A0_coriolisvec_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-20 17:54:29
% EndTime: 2018-12-20 17:54:30
% DurationCPUTime: 1.00s
% Computational Cost: add. (3277->166), mult. (4877->312), div. (666->3), fcn. (3892->14), ass. (0->149)
t406 = 2 * mrSges(2,3);
t355 = xDP(3);
t344 = t355 ^ 2;
t349 = sin(qJ(1,3));
t352 = cos(qJ(1,3));
t358 = pkin(1) + pkin(2);
t325 = -t352 * qJ(2,3) + t349 * t358;
t328 = t349 * qJ(2,3) + t358 * t352;
t346 = legFrame(3,3);
t334 = sin(t346);
t337 = cos(t346);
t298 = t325 * t337 + t334 * t328;
t301 = -t334 * t325 + t328 * t337;
t359 = xP(3);
t341 = sin(t359);
t342 = cos(t359);
t368 = koppelP(3,2);
t371 = koppelP(3,1);
t322 = -t341 * t368 + t342 * t371;
t356 = xDP(2);
t304 = t322 * t355 + t356;
t319 = t341 * t371 + t342 * t368;
t357 = xDP(1);
t307 = -t319 * t355 + t357;
t363 = 0.1e1 / qJ(2,3);
t265 = (t298 * t304 + t301 * t307) * t363;
t405 = 0.2e1 * t265;
t350 = sin(qJ(1,2));
t353 = cos(qJ(1,2));
t326 = -t353 * qJ(2,2) + t350 * t358;
t329 = t350 * qJ(2,2) + t358 * t353;
t347 = legFrame(2,3);
t335 = sin(t347);
t338 = cos(t347);
t299 = t326 * t338 + t335 * t329;
t302 = -t335 * t326 + t329 * t338;
t369 = koppelP(2,2);
t372 = koppelP(2,1);
t323 = -t341 * t369 + t342 * t372;
t305 = t323 * t355 + t356;
t320 = t341 * t372 + t342 * t369;
t308 = -t320 * t355 + t357;
t365 = 0.1e1 / qJ(2,2);
t266 = (t299 * t305 + t302 * t308) * t365;
t404 = 0.2e1 * t266;
t351 = sin(qJ(1,1));
t354 = cos(qJ(1,1));
t327 = -t354 * qJ(2,1) + t351 * t358;
t330 = t351 * qJ(2,1) + t358 * t354;
t348 = legFrame(1,3);
t336 = sin(t348);
t339 = cos(t348);
t300 = t327 * t339 + t336 * t330;
t303 = -t336 * t327 + t330 * t339;
t370 = koppelP(1,2);
t373 = koppelP(1,1);
t324 = -t341 * t370 + t342 * t373;
t306 = t324 * t355 + t356;
t321 = t341 * t373 + t342 * t370;
t309 = -t321 * t355 + t357;
t367 = 0.1e1 / qJ(2,1);
t267 = (t300 * t306 + t303 * t309) * t367;
t403 = 0.2e1 * t267;
t402 = 0.2e1 * t358;
t310 = t334 * t352 + t337 * t349;
t311 = -t334 * t349 + t337 * t352;
t274 = (t304 * t310 + t307 * t311) * t363;
t331 = m(2) * qJ(2,3) + mrSges(2,3);
t401 = t274 ^ 2 * t331;
t312 = t335 * t353 + t338 * t350;
t313 = -t335 * t350 + t338 * t353;
t275 = (t305 * t312 + t308 * t313) * t365;
t332 = m(2) * qJ(2,2) + mrSges(2,3);
t400 = t275 ^ 2 * t332;
t314 = t336 * t354 + t339 * t351;
t315 = -t336 * t351 + t339 * t354;
t276 = (t306 * t314 + t309 * t315) * t367;
t333 = m(2) * qJ(2,1) + mrSges(2,3);
t399 = t276 ^ 2 * t333;
t398 = t274 * t363;
t397 = t275 * t365;
t396 = t276 * t367;
t362 = qJ(2,3) ^ 2;
t375 = (pkin(1) ^ 2);
t383 = -t375 + (-2 * pkin(1) - pkin(2)) * pkin(2);
t253 = (t265 * t402 + (-t362 + t383) * t274) * t398;
t256 = (-t358 * t274 + t405) * t398;
t379 = 2 * pkin(1) * mrSges(2,1) + Ifges(2,2) + Ifges(1,3);
t316 = m(2) * (t362 + t375) + qJ(2,3) * t406 + t379;
t340 = m(2) * pkin(1) + mrSges(2,1);
t247 = t340 * t253 - t316 * t256;
t395 = t363 * t247;
t250 = -m(2) * t253 + t340 * t256;
t394 = t363 * t250;
t393 = t363 * t344;
t364 = qJ(2,2) ^ 2;
t254 = (t266 * t402 + (-t364 + t383) * t275) * t397;
t257 = (-t358 * t275 + t404) * t397;
t317 = m(2) * (t364 + t375) + qJ(2,2) * t406 + t379;
t248 = t340 * t254 - t317 * t257;
t392 = t365 * t248;
t251 = -m(2) * t254 + t340 * t257;
t391 = t365 * t251;
t390 = t365 * t344;
t366 = qJ(2,1) ^ 2;
t255 = (t267 * t402 + (-t366 + t383) * t276) * t396;
t258 = (-t358 * t276 + t403) * t396;
t318 = m(2) * (t366 + t375) + qJ(2,1) * t406 + t379;
t249 = t340 * t255 - t318 * t258;
t389 = t367 * t249;
t252 = -m(2) * t255 + t340 * t258;
t388 = t367 * t252;
t387 = t367 * t344;
t386 = t363 * t401;
t385 = t365 * t400;
t384 = t367 * t399;
t382 = t274 * t331 * t405;
t381 = t275 * t332 * t404;
t380 = t276 * t333 * t403;
t378 = t363 * t382;
t377 = t365 * t381;
t376 = t367 * t380;
t361 = mrSges(3,1);
t360 = mrSges(3,2);
t291 = (m(2) * t303 - t315 * t340) * t367;
t290 = (m(2) * t300 - t314 * t340) * t367;
t289 = (m(2) * t302 - t313 * t340) * t365;
t288 = (m(2) * t299 - t312 * t340) * t365;
t287 = (m(2) * t301 - t311 * t340) * t363;
t286 = (m(2) * t298 - t310 * t340) * t363;
t285 = (t314 * t324 - t315 * t321) * t367;
t284 = (t312 * t323 - t313 * t320) * t365;
t283 = (t310 * t322 - t311 * t319) * t363;
t282 = (-t303 * t340 + t315 * t318) * t367;
t281 = (-t300 * t340 + t314 * t318) * t367;
t280 = (-t302 * t340 + t313 * t317) * t365;
t279 = (-t299 * t340 + t312 * t317) * t365;
t278 = (-t301 * t340 + t311 * t316) * t363;
t277 = (-t298 * t340 + t310 * t316) * t363;
t270 = (t300 * t324 - t303 * t321) * t367;
t269 = (t299 * t323 - t302 * t320) * t365;
t268 = (t298 * t322 - t301 * t319) * t363;
t264 = t270 * m(2) - t285 * t340;
t263 = t269 * m(2) - t284 * t340;
t262 = t268 * m(2) - t283 * t340;
t261 = -t270 * t340 + t285 * t318;
t260 = -t269 * t340 + t284 * t317;
t259 = -t268 * t340 + t283 * t316;
t1 = [(-(t282 * t315 + t291 * t303) * t324 - (t282 * t314 + t291 * t300) * t321) * t387 + t315 * t389 + t303 * t388 + t315 * t376 - t303 * t384 + (-(t280 * t313 + t289 * t302) * t323 - (t280 * t312 + t289 * t299) * t320) * t390 + t313 * t392 + t302 * t391 + t313 * t377 - t302 * t385 + (-(t278 * t311 + t287 * t301) * t322 - (t278 * t310 + t287 * t298) * t319) * t393 + t311 * t395 + t301 * t394 + t311 * t378 - t301 * t386 - t344 * (-t341 * t360 + t342 * t361); (-(t281 * t315 + t290 * t303) * t324 - (t281 * t314 + t290 * t300) * t321) * t387 + t314 * t389 + t300 * t388 + t314 * t376 - t300 * t384 + (-(t279 * t313 + t288 * t302) * t323 - (t279 * t312 + t288 * t299) * t320) * t390 + t312 * t392 + t299 * t391 + t312 * t377 - t299 * t385 + (-(t277 * t311 + t286 * t301) * t322 - (t277 * t310 + t286 * t298) * t319) * t393 + t310 * t395 + t298 * t394 + t310 * t378 - t298 * t386 - t344 * (t341 * t361 + t342 * t360); (-(t261 * t315 + t264 * t303) * t324 - (t261 * t314 + t264 * t300) * t321) * t387 + (-(t260 * t313 + t263 * t302) * t323 - (t260 * t312 + t263 * t299) * t320) * t390 + (-(t259 * t311 + t262 * t301) * t322 - (t259 * t310 + t262 * t298) * t319) * t393 + (t249 + t380) * t285 + (t248 + t381) * t284 + (t247 + t382) * t283 + (t252 - t399) * t270 + (t251 - t400) * t269 + (t250 - t401) * t268;];
taucX  = t1;
