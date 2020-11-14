% Calculate inertia matrix for parallel robot
% P3PRRRR1G1A0
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% Datum: 2020-03-09 20:34
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3PRRRR1G1A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR1G1A0_inertia_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR1G1A0_inertia_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRRR1G1A0_inertia_para_pf_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR1G1A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRRR1G1A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3PRRRR1G1A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR1G1A0_inertia_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR1G1A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 20:34:02
% EndTime: 2020-03-09 20:34:02
% DurationCPUTime: 0.62s
% Computational Cost: add. (741->128), mult. (1563->262), div. (408->10), fcn. (1260->24), ass. (0->158)
t446 = m(3) / 0.2e1;
t445 = Icges(3,2) / 0.2e1;
t444 = rSges(3,3) * m(3);
t370 = rSges(3,2) ^ 2;
t371 = rSges(3,1) ^ 2;
t443 = (-t370 + t371) * t446 - Icges(3,1) / 0.2e1 + t445;
t352 = sin(qJ(3,3));
t358 = cos(qJ(3,3));
t322 = t352 * rSges(3,1) + t358 * rSges(3,2);
t442 = m(3) * t322;
t354 = sin(qJ(3,2));
t360 = cos(qJ(3,2));
t323 = t354 * rSges(3,1) + t360 * rSges(3,2);
t441 = m(3) * t323;
t356 = sin(qJ(3,1));
t362 = cos(qJ(3,1));
t324 = t356 * rSges(3,1) + t362 * rSges(3,2);
t440 = m(3) * t324;
t372 = 0.1e1 / pkin(2);
t439 = m(3) * t372;
t329 = m(2) * rSges(2,2) - t444;
t353 = sin(qJ(2,3));
t359 = cos(qJ(2,3));
t366 = m(2) * rSges(2,1);
t301 = (t366 + (rSges(3,1) * t358 - rSges(3,2) * t352) * m(3)) * t359 - t353 * t329;
t343 = 0.1e1 / t358;
t438 = t301 * t343;
t355 = sin(qJ(2,2));
t361 = cos(qJ(2,2));
t302 = (t366 + (rSges(3,1) * t360 - rSges(3,2) * t354) * m(3)) * t361 - t355 * t329;
t345 = 0.1e1 / t360;
t437 = t302 * t345;
t357 = sin(qJ(2,1));
t363 = cos(qJ(2,1));
t303 = (t366 + (rSges(3,1) * t362 - rSges(3,2) * t356) * m(3)) * t363 - t357 * t329;
t347 = 0.1e1 / t362;
t436 = t303 * t347;
t410 = t370 + t371;
t435 = (t410 * m(3) + Icges(3,3)) * t372;
t349 = legFrame(3,3);
t333 = sin(t349);
t434 = t333 * t343;
t350 = legFrame(2,3);
t334 = sin(t350);
t433 = t334 * t345;
t351 = legFrame(1,3);
t335 = sin(t351);
t432 = t335 * t347;
t336 = cos(t349);
t431 = t336 * t343;
t337 = cos(t350);
t430 = t337 * t345;
t338 = cos(t351);
t429 = t338 * t347;
t340 = 0.1e1 / t353;
t428 = t340 * t343;
t344 = 0.1e1 / t358 ^ 2;
t427 = t340 * t344;
t341 = 0.1e1 / t355;
t426 = t341 * t345;
t346 = 0.1e1 / t360 ^ 2;
t425 = t341 * t346;
t342 = 0.1e1 / t357;
t424 = t342 * t347;
t348 = 0.1e1 / t362 ^ 2;
t423 = t342 * t348;
t422 = t344 * t372;
t421 = t346 * t372;
t420 = t348 * t372;
t419 = t352 * t359;
t418 = t354 * t361;
t417 = t356 * t363;
t416 = t358 * t359;
t415 = t360 * t361;
t414 = t362 * t363;
t307 = -t333 * t419 - t336 * t358;
t310 = t333 * t352 + t336 * t416;
t385 = t340 * t422;
t382 = t301 * t385;
t339 = m(1) + m(2) + m(3);
t388 = t339 * t428;
t413 = t307 * t382 + t310 * t388;
t308 = -t334 * t418 - t337 * t360;
t312 = t334 * t354 + t337 * t415;
t384 = t341 * t421;
t381 = t302 * t384;
t387 = t339 * t426;
t412 = t308 * t381 + t312 * t387;
t309 = -t335 * t417 - t338 * t362;
t314 = t335 * t356 + t338 * t414;
t383 = t342 * t420;
t380 = t303 * t383;
t386 = t339 * t424;
t411 = t309 * t380 + t314 * t386;
t332 = -m(3) * rSges(3,1) * rSges(3,2) + Icges(3,4);
t367 = 0.2e1 * qJ(3,3);
t373 = Icges(2,3) + (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + (0.2e1 * rSges(3,3) ^ 2 + t410) * t446 + t445 + Icges(3,1) / 0.2e1;
t409 = (cos(t367) * t443 + t332 * sin(t367) + t373) * t422;
t368 = 0.2e1 * qJ(3,2);
t408 = (cos(t368) * t443 + t332 * sin(t368) + t373) * t421;
t369 = 0.2e1 * qJ(3,1);
t407 = (cos(t369) * t443 + t332 * sin(t369) + t373) * t420;
t406 = t307 * t427;
t405 = t308 * t425;
t404 = t309 * t423;
t403 = t310 * t428;
t311 = -t358 * t333 + t336 * t419;
t402 = t311 * t427;
t401 = t312 * t426;
t313 = -t360 * t334 + t337 * t418;
t400 = t313 * t425;
t399 = t314 * t424;
t315 = -t362 * t335 + t338 * t417;
t398 = t315 * t423;
t316 = t333 * t416 - t336 * t352;
t397 = t316 * t428;
t317 = t334 * t415 - t337 * t354;
t396 = t317 * t426;
t318 = t335 * t414 - t338 * t356;
t395 = t318 * t424;
t330 = -rSges(3,2) * t444 + Icges(3,6);
t331 = rSges(3,1) * t444 - Icges(3,5);
t319 = t330 * t358 - t331 * t352;
t394 = t319 * t343 * t372;
t320 = t330 * t360 - t331 * t354;
t393 = t320 * t345 * t372;
t321 = t330 * t362 - t331 * t356;
t392 = t321 * t347 * t372;
t391 = t322 * t343 * t353;
t390 = t323 * t345 * t355;
t389 = t324 * t347 * t357;
t376 = t391 * t439;
t277 = t311 * t382 + t316 * t388 + t336 * t376;
t375 = t390 * t439;
t278 = t313 * t381 + t317 * t387 + t337 * t375;
t374 = t389 * t439;
t279 = t315 * t380 + t318 * t386 + t338 * t374;
t379 = t319 * t385;
t378 = t320 * t384;
t377 = t321 * t383;
t285 = t315 * t377 + (-t318 * t440 - t338 * t435) * t347;
t284 = t313 * t378 + (-t317 * t441 - t337 * t435) * t345;
t283 = t311 * t379 + (-t316 * t442 - t336 * t435) * t343;
t282 = t309 * t377 + (-t314 * t440 + t335 * t435) * t347;
t281 = t308 * t378 + (-t312 * t441 + t334 * t435) * t345;
t280 = t307 * t379 + (-t310 * t442 + t333 * t435) * t343;
t276 = -t335 * t374 + t411;
t275 = -t334 * t375 + t412;
t274 = -t333 * t376 + t413;
t273 = -t338 * t392 + (t315 * t407 + t318 * t436) * t342;
t272 = -t337 * t393 + (t313 * t408 + t317 * t437) * t341;
t271 = -t336 * t394 + (t311 * t409 + t316 * t438) * t340;
t270 = t335 * t392 + (t309 * t407 + t314 * t436) * t342;
t269 = t334 * t393 + (t308 * t408 + t312 * t437) * t341;
t268 = t333 * t394 + (t307 * t409 + t310 * t438) * t340;
t267 = t279 + t278 + t277;
t266 = (-t333 * t391 - t334 * t390 - t335 * t389) * t439 + t411 + t412 + t413;
t1 = [t274 * t403 + t275 * t401 + t276 * t399 + m(4) + (t268 * t406 + t269 * t405 + t270 * t404 + t280 * t434 + t281 * t433 + t282 * t432) * t372, t274 * t397 + t275 * t396 + t276 * t395 + (t268 * t402 + t269 * t400 + t270 * t398 - t280 * t431 - t281 * t430 - t282 * t429) * t372, t266; t277 * t403 + t278 * t401 + t279 * t399 + (t271 * t406 + t272 * t405 + t273 * t404 + t283 * t434 + t284 * t433 + t285 * t432) * t372, t277 * t397 + t278 * t396 + t279 * t395 + m(4) + (t271 * t402 + t272 * t400 + t273 * t398 - t283 * t431 - t284 * t430 - t285 * t429) * t372, t267; t266, t267, 0.3e1 * m(1) + 0.3e1 * m(2) + 0.3e1 * m(3) + m(4);];
MX  = t1;
