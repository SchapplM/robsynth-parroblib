% Calculate inertia matrix for parallel robot
% P3RPRRR12V1G2A0
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% Datum: 2020-08-06 18:25
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RPRRR12V1G2A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR12V1G2A0_inertia_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR12V1G2A0_inertia_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RPRRR12V1G2A0_inertia_para_pf_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR12V1G2A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RPRRR12V1G2A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RPRRR12V1G2A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR12V1G2A0_inertia_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR12V1G2A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:24:27
% EndTime: 2020-08-06 18:24:28
% DurationCPUTime: 0.92s
% Computational Cost: add. (2241->204), mult. (2913->382), div. (396->7), fcn. (1671->18), ass. (0->154)
t454 = 2 * rSges(2,3);
t453 = 2 * m(3) * rSges(3,1) * rSges(3,2) - 2 * Icges(3,4);
t384 = (rSges(3,3) + pkin(5));
t372 = sin(qJ(3,3));
t378 = cos(qJ(3,3));
t452 = m(3) * (t378 * rSges(3,1) - t372 * rSges(3,2));
t374 = sin(qJ(3,2));
t380 = cos(qJ(3,2));
t451 = m(3) * (t380 * rSges(3,1) - t374 * rSges(3,2));
t376 = sin(qJ(3,1));
t382 = cos(qJ(3,1));
t450 = m(3) * (t382 * rSges(3,1) - t376 * rSges(3,2));
t449 = pkin(3) * t378;
t379 = cos(qJ(1,3));
t448 = pkin(3) * t379;
t447 = pkin(3) * t380;
t381 = cos(qJ(1,2));
t446 = pkin(3) * t381;
t445 = pkin(3) * t382;
t383 = cos(qJ(1,1));
t444 = pkin(3) * t383;
t443 = (-pkin(1) - t384) * m(3);
t360 = pkin(1) + pkin(5) + pkin(6);
t373 = sin(qJ(1,3));
t333 = qJ(2,3) * t379 - t360 * t373;
t369 = legFrame(3,2);
t353 = sin(t369);
t356 = cos(t369);
t412 = t378 * qJ(2,3);
t318 = (-t333 * t356 + t353 * t449) * t372 + (t378 - 0.1e1) * (t378 + 0.1e1) * t356 * t448 + t353 * t412;
t362 = 0.1e1 / t372;
t442 = t318 * t362;
t375 = sin(qJ(1,2));
t334 = qJ(2,2) * t381 - t360 * t375;
t370 = legFrame(2,2);
t354 = sin(t370);
t357 = cos(t370);
t411 = t380 * qJ(2,2);
t319 = (-t334 * t357 + t354 * t447) * t374 + (t380 - 0.1e1) * (t380 + 0.1e1) * t357 * t446 + t354 * t411;
t363 = 0.1e1 / t374;
t441 = t319 * t363;
t377 = sin(qJ(1,1));
t335 = qJ(2,1) * t383 - t360 * t377;
t371 = legFrame(1,2);
t355 = sin(t371);
t358 = cos(t371);
t410 = t382 * qJ(2,1);
t320 = (-t335 * t358 + t355 * t445) * t376 + (t382 - 0.1e1) * (t382 + 0.1e1) * t358 * t444 + t355 * t410;
t364 = 0.1e1 / t376;
t440 = t320 * t364;
t365 = t378 ^ 2;
t321 = (t333 * t353 + t356 * t449) * t372 + (-t365 + 0.1e1) * t353 * t448 + t356 * t412;
t439 = t321 * t362;
t366 = t380 ^ 2;
t322 = (t334 * t354 + t357 * t447) * t374 + (-t366 + 0.1e1) * t354 * t446 + t357 * t411;
t438 = t322 * t363;
t367 = t382 ^ 2;
t323 = (t335 * t355 + t358 * t445) * t376 + (-t367 + 0.1e1) * t355 * t444 + t358 * t410;
t437 = t323 * t364;
t347 = t372 * pkin(3) + qJ(2,3);
t330 = t373 * t347 + t360 * t379;
t336 = t443 - (pkin(1) - rSges(2,2)) * m(2);
t344 = 0.1e1 / t347;
t385 = m(2) + m(3);
t324 = (t330 * t385 + t336 * t379) * t344;
t436 = t324 * t362;
t348 = t374 * pkin(3) + qJ(2,2);
t331 = t375 * t348 + t360 * t381;
t345 = 0.1e1 / t348;
t325 = (t331 * t385 + t336 * t381) * t345;
t435 = t325 * t363;
t349 = t376 * pkin(3) + qJ(2,1);
t332 = t377 * t349 + t360 * t383;
t346 = 0.1e1 / t349;
t326 = (t332 * t385 + t336 * t383) * t346;
t434 = t326 * t364;
t433 = t336 * t362;
t432 = t336 * t363;
t431 = t336 * t364;
t430 = t353 * t362;
t429 = t353 * t373;
t428 = t354 * t363;
t427 = t354 * t375;
t426 = t355 * t364;
t425 = t355 * t377;
t424 = t356 * t362;
t423 = t356 * t373;
t422 = t357 * t363;
t421 = t357 * t375;
t420 = t358 * t364;
t419 = t358 * t377;
t418 = t362 * t385;
t393 = 0.1e1 / pkin(3);
t417 = t362 * t393;
t416 = t363 * t385;
t415 = t363 * t393;
t414 = t364 * t385;
t413 = t364 * t393;
t409 = t362 * t452;
t408 = t363 * t451;
t407 = t364 * t450;
t406 = t353 * t417;
t405 = t354 * t415;
t404 = t355 * t413;
t403 = t356 * t417;
t402 = t357 * t415;
t401 = t358 * t413;
t400 = t393 * t409;
t399 = t393 * t408;
t398 = t393 * t407;
t397 = (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + Icges(2,1) + Icges(3,2) + Icges(1,3);
t392 = rSges(3,1) ^ 2;
t394 = pkin(1) ^ 2;
t396 = t392 + t394 + (2 * pkin(1) + t384) * t384;
t395 = rSges(2,3) ^ 2 + t394 + (-2 * pkin(1) + rSges(2,2)) * rSges(2,2);
t390 = rSges(3,2) ^ 2;
t388 = qJ(2,1) ^ 2;
t387 = qJ(2,2) ^ 2;
t386 = qJ(2,3) ^ 2;
t343 = (t390 + t392) * m(3) + Icges(3,3);
t342 = rSges(3,1) * t443 + Icges(3,5);
t341 = -rSges(3,2) * t443 - Icges(3,6);
t340 = (t390 - t392) * m(3) + Icges(3,1) - Icges(3,2);
t329 = t341 * t376 + t342 * t382;
t328 = t341 * t374 + t342 * t380;
t327 = t341 * t372 + t342 * t378;
t317 = (t329 * t383 + t332 * t450) * t346;
t316 = (t328 * t381 + t331 * t451) * t345;
t315 = (t327 * t379 + t330 * t452) * t344;
t314 = t340 * t367 + t376 * t382 * t453 + (qJ(2,1) * t454 + t388 + t395) * m(2) + (t388 + 0.2e1 * (rSges(3,1) * t376 + rSges(3,2) * t382) * qJ(2,1) + t396) * m(3) + t397;
t313 = t340 * t366 + t374 * t380 * t453 + (qJ(2,2) * t454 + t387 + t395) * m(2) + (t387 + 0.2e1 * (rSges(3,1) * t374 + rSges(3,2) * t380) * qJ(2,2) + t396) * m(3) + t397;
t312 = t340 * t365 + t372 * t378 * t453 + (qJ(2,3) * t454 + t386 + t395) * m(2) + (t386 + 0.2e1 * (rSges(3,1) * t372 + rSges(3,2) * t378) * qJ(2,3) + t396) * m(3) + t397;
t311 = -t358 * t398 + (t323 * t414 - t336 * t425) * t346;
t310 = -t357 * t399 + (t322 * t416 - t336 * t427) * t345;
t309 = -t356 * t400 + (t321 * t418 - t336 * t429) * t344;
t308 = -t355 * t398 + (t320 * t414 + t336 * t419) * t346;
t307 = -t354 * t399 + (t319 * t416 + t336 * t421) * t345;
t306 = -t353 * t400 + (t318 * t418 + t336 * t423) * t344;
t305 = (t314 * t383 + t332 * t336) * t346;
t304 = (t313 * t381 + t331 * t336) * t345;
t303 = (t312 * t379 + t330 * t336) * t344;
t302 = -t343 * t401 + (t323 * t407 - t329 * t425) * t346;
t301 = -t343 * t402 + (t322 * t408 - t328 * t427) * t345;
t300 = -t343 * t403 + (t321 * t409 - t327 * t429) * t344;
t299 = -t343 * t404 + (t320 * t407 + t329 * t419) * t346;
t298 = -t343 * t405 + (t319 * t408 + t328 * t421) * t345;
t297 = -t343 * t406 + (t318 * t409 + t327 * t423) * t344;
t296 = -t329 * t401 + (-t314 * t425 + t323 * t431) * t346;
t295 = -t328 * t402 + (-t313 * t427 + t322 * t432) * t345;
t294 = -t327 * t403 + (-t312 * t429 + t321 * t433) * t344;
t293 = -t329 * t404 + (t314 * t419 + t320 * t431) * t346;
t292 = -t328 * t405 + (t313 * t421 + t319 * t432) * t345;
t291 = -t327 * t406 + (t312 * t423 + t318 * t433) * t344;
t1 = [m(4) + (t293 * t419 + t308 * t440) * t346 + (t292 * t421 + t307 * t441) * t345 + (t291 * t423 + t306 * t442) * t344 + (-t297 * t430 - t298 * t428 - t299 * t426) * t393, (-t293 * t425 + t308 * t437) * t346 + (-t292 * t427 + t307 * t438) * t345 + (-t291 * t429 + t306 * t439) * t344 + (-t297 * t424 - t298 * t422 - t299 * t420) * t393, (t293 * t383 + t308 * t332) * t346 + (t292 * t381 + t307 * t331) * t345 + (t291 * t379 + t306 * t330) * t344; (t296 * t419 + t311 * t440) * t346 + (t295 * t421 + t310 * t441) * t345 + (t294 * t423 + t309 * t442) * t344 + (-t300 * t430 - t301 * t428 - t302 * t426) * t393, m(4) + (-t296 * t425 + t311 * t437) * t346 + (-t295 * t427 + t310 * t438) * t345 + (-t294 * t429 + t309 * t439) * t344 + (-t300 * t424 - t301 * t422 - t302 * t420) * t393, (t296 * t383 + t311 * t332) * t346 + (t295 * t381 + t310 * t331) * t345 + (t294 * t379 + t309 * t330) * t344; (t305 * t419 + t320 * t434) * t346 + (t304 * t421 + t319 * t435) * t345 + (t303 * t423 + t318 * t436) * t344 + (-t315 * t430 - t316 * t428 - t317 * t426) * t393, (-t305 * t425 + t323 * t434) * t346 + (-t304 * t427 + t322 * t435) * t345 + (-t303 * t429 + t321 * t436) * t344 + (-t315 * t424 - t316 * t422 - t317 * t420) * t393, m(4) + (t305 * t383 + t326 * t332) * t346 + (t304 * t381 + t325 * t331) * t345 + (t303 * t379 + t324 * t330) * t344;];
MX  = t1;
