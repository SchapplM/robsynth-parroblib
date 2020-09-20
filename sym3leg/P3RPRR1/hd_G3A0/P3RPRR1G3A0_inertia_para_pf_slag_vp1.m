% Calculate inertia matrix for parallel robot
% P3RPRR1G3P3A0
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
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Icges [4x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
%
% Output:
% MX [3x3]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:27
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RPRR1G3P3A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRR1G3P3A0_inertia_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRR1G3P3A0_inertia_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRR1G3P3A0_inertia_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRR1G3P3A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RPRR1G3P3A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RPRR1G3P3A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRR1G3P3A0_inertia_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRR1G3P3A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:26:46
% EndTime: 2020-03-09 21:26:46
% DurationCPUTime: 0.75s
% Computational Cost: add. (2961->170), mult. (3048->225), div. (243->4), fcn. (1260->65), ass. (0->118)
t366 = pkin(7) + qJ(3,3);
t394 = pkin(1) * sin(t366) + sin(qJ(3,3)) * pkin(2);
t308 = 0.1e1 / t394;
t405 = t308 / 0.2e1;
t367 = pkin(7) + qJ(3,2);
t393 = pkin(1) * sin(t367) + sin(qJ(3,2)) * pkin(2);
t309 = 0.1e1 / t393;
t403 = t309 / 0.2e1;
t368 = pkin(7) + qJ(3,1);
t392 = pkin(1) * sin(t368) + sin(qJ(3,1)) * pkin(2);
t310 = 0.1e1 / t392;
t401 = t310 / 0.2e1;
t383 = 0.1e1 / pkin(3);
t434 = t383 / 0.2e1;
t371 = legFrame(1,2);
t398 = qJ(1,1) + pkin(7);
t340 = t371 + t398;
t333 = qJ(3,1) + t340;
t341 = -t371 + t398;
t334 = qJ(3,1) + t341;
t307 = cos(t334) + cos(t333);
t433 = t307 * t401;
t370 = legFrame(2,2);
t397 = qJ(1,2) + pkin(7);
t338 = t370 + t397;
t331 = qJ(3,2) + t338;
t339 = -t370 + t397;
t332 = qJ(3,2) + t339;
t306 = cos(t332) + cos(t331);
t432 = t306 * t403;
t369 = legFrame(3,2);
t396 = qJ(1,3) + pkin(7);
t336 = t369 + t396;
t329 = qJ(3,3) + t336;
t337 = -t369 + t396;
t330 = qJ(3,3) + t337;
t305 = cos(t330) + cos(t329);
t431 = t305 * t405;
t304 = -sin(t333) + sin(t334);
t430 = t304 * t401;
t303 = -sin(t331) + sin(t332);
t429 = t303 * t403;
t302 = -sin(t329) + sin(t330);
t428 = t302 * t405;
t427 = (-t392 * rSges(3,2) + (pkin(1) * cos(t368) + pkin(2) * cos(qJ(3,1))) * rSges(3,1)) * m(3);
t426 = (-t393 * rSges(3,2) + (pkin(1) * cos(t367) + cos(qJ(3,2)) * pkin(2)) * rSges(3,1)) * m(3);
t425 = (-t394 * rSges(3,2) + (pkin(1) * cos(t366) + cos(qJ(3,3)) * pkin(2)) * rSges(3,1)) * m(3);
t354 = qJ(1,3) + t369;
t355 = qJ(1,3) - t369;
t290 = t302 * pkin(3) + (-sin(t336) + sin(t337)) * pkin(2) + (-sin(t354) + sin(t355)) * pkin(1);
t424 = t290 * t308;
t356 = qJ(1,2) + t370;
t357 = qJ(1,2) - t370;
t291 = t303 * pkin(3) + (-sin(t338) + sin(t339)) * pkin(2) + (-sin(t356) + sin(t357)) * pkin(1);
t423 = t291 * t309;
t358 = qJ(1,1) + t371;
t359 = qJ(1,1) - t371;
t292 = t304 * pkin(3) + (sin(t341) - sin(t340)) * pkin(2) + (-sin(t358) + sin(t359)) * pkin(1);
t422 = t292 * t310;
t293 = -t305 * pkin(3) + (-cos(t337) - cos(t336)) * pkin(2) + (-cos(t355) - cos(t354)) * pkin(1);
t421 = t293 * t308;
t294 = -t306 * pkin(3) + (-cos(t339) - cos(t338)) * pkin(2) + (-cos(t357) - cos(t356)) * pkin(1);
t420 = t294 * t309;
t295 = -t307 * pkin(3) + (-cos(t341) - cos(t340)) * pkin(2) + (-cos(t359) - cos(t358)) * pkin(1);
t419 = t295 * t310;
t399 = rSges(3,1) ^ 2 + rSges(3,2) ^ 2;
t328 = t399 * m(3) + Icges(3,3);
t296 = t425 + t328;
t418 = t296 * t383;
t297 = t426 + t328;
t417 = t297 * t383;
t298 = t427 + t328;
t416 = t298 * t383;
t342 = sin(qJ(1,3) + t366);
t299 = pkin(3) * t342 + pkin(2) * sin(t396) + sin(qJ(1,3)) * pkin(1);
t415 = t299 * t308;
t343 = sin(qJ(1,2) + t367);
t300 = pkin(3) * t343 + pkin(2) * sin(t397) + sin(qJ(1,2)) * pkin(1);
t414 = t300 * t309;
t344 = sin(qJ(1,1) + t368);
t301 = pkin(3) * t344 + pkin(2) * sin(t398) + sin(qJ(1,1)) * pkin(1);
t413 = t301 * t310;
t406 = t308 * t342;
t404 = t309 * t343;
t402 = t310 * t344;
t400 = t328 * t383;
t360 = sin(t369);
t361 = sin(t370);
t362 = sin(t371);
t363 = cos(t369);
t364 = cos(t370);
t365 = cos(t371);
t379 = m(2) + m(3);
t395 = (t360 * t363 + t361 * t364 + t362 * t365) * t379;
t384 = pkin(1) ^ 2;
t385 = Icges(1,3) + Icges(2,3) + Icges(3,3) + (pkin(2) ^ 2 + t384 + t399) * m(3) + (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + (rSges(2,1) ^ 2 + rSges(2,2) ^ 2 + t384) * m(2) + 0.2e1 * ((m(2) * rSges(2,1) + m(3) * pkin(2)) * cos(pkin(7)) - m(2) * sin(pkin(7)) * rSges(2,2)) * pkin(1);
t289 = t385 + 0.2e1 * t427;
t288 = t385 + 0.2e1 * t426;
t287 = t385 + 0.2e1 * t425;
t286 = (-t298 * t344 + t301 * t400) * t310;
t285 = (-t297 * t343 + t300 * t400) * t309;
t284 = (-t296 * t342 + t299 * t400) * t308;
t283 = (t295 * t400 + t298 * t307) * t401;
t282 = (t294 * t400 + t297 * t306) * t403;
t281 = (t293 * t400 + t296 * t305) * t405;
t280 = (-t292 * t400 + t298 * t304) * t401;
t279 = (-t291 * t400 + t297 * t303) * t403;
t278 = (-t290 * t400 + t296 * t302) * t405;
t277 = (-t289 * t344 + t301 * t416) * t310;
t276 = (-t288 * t343 + t300 * t417) * t309;
t275 = (-t287 * t342 + t299 * t418) * t308;
t274 = (t289 * t307 + t295 * t416) * t401;
t273 = (t288 * t306 + t294 * t417) * t403;
t272 = (t287 * t305 + t293 * t418) * t405;
t271 = (t289 * t304 - t292 * t416) * t401;
t270 = (t288 * t303 - t291 * t417) * t403;
t269 = (t287 * t302 - t290 * t418) * t405;
t1 = [m(4) + (t360 ^ 2 + t361 ^ 2 + t362 ^ 2) * t379 + t272 * t431 + t273 * t432 + t274 * t433 + (t281 * t421 + t282 * t420 + t283 * t419) * t434, t272 * t428 + t273 * t429 + t274 * t430 + (-t281 * t424 - t282 * t423 - t283 * t422) * t434 + t395, -t272 * t406 - t273 * t404 - t274 * t402 + (t281 * t415 + t282 * t414 + t283 * t413) * t383; t269 * t431 + t270 * t432 + t271 * t433 + (t278 * t421 + t279 * t420 + t280 * t419) * t434 + t395, m(4) + (t363 ^ 2 + t364 ^ 2 + t365 ^ 2) * t379 + t269 * t428 + t270 * t429 + t271 * t430 + (-t278 * t424 - t279 * t423 - t280 * t422) * t434, -t269 * t406 - t270 * t404 - t271 * t402 + (t278 * t415 + t279 * t414 + t280 * t413) * t383; t275 * t431 + t276 * t432 + t277 * t433 + (t284 * t421 + t285 * t420 + t286 * t419) * t434, t275 * t428 + t276 * t429 + t277 * t430 + (-t284 * t424 - t285 * t423 - t286 * t422) * t434, -t275 * t406 - t276 * t404 - t277 * t402 + m(4) + (t284 * t415 + t285 * t414 + t286 * t413) * t383;];
MX  = t1;
