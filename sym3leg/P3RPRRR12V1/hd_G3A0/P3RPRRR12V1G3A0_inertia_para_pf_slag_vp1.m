% Calculate inertia matrix for parallel robot
% P3RPRRR12V1G3A0
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
% Datum: 2020-08-06 18:28
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RPRRR12V1G3A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR12V1G3A0_inertia_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR12V1G3A0_inertia_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RPRRR12V1G3A0_inertia_para_pf_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR12V1G3A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RPRRR12V1G3A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RPRRR12V1G3A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR12V1G3A0_inertia_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR12V1G3A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:27:59
% EndTime: 2020-08-06 18:28:00
% DurationCPUTime: 1.05s
% Computational Cost: add. (2241->204), mult. (2946->385), div. (396->7), fcn. (1704->18), ass. (0->151)
t447 = 2 * rSges(2,3);
t446 = 2 * m(3) * rSges(3,1) * rSges(3,2) - 2 * Icges(3,4);
t377 = (rSges(3,3) + pkin(5));
t365 = sin(qJ(3,3));
t371 = cos(qJ(3,3));
t445 = m(3) * (t371 * rSges(3,1) - t365 * rSges(3,2));
t367 = sin(qJ(3,2));
t373 = cos(qJ(3,2));
t444 = m(3) * (t373 * rSges(3,1) - t367 * rSges(3,2));
t369 = sin(qJ(3,1));
t375 = cos(qJ(3,1));
t443 = m(3) * (t375 * rSges(3,1) - t369 * rSges(3,2));
t442 = (-pkin(1) - t377) * m(3);
t362 = legFrame(3,2);
t346 = sin(t362);
t349 = cos(t362);
t358 = t371 ^ 2;
t366 = sin(qJ(1,3));
t353 = pkin(1) + pkin(5) + pkin(6);
t372 = cos(qJ(1,3));
t388 = qJ(2,3) * t366 + t353 * t372;
t408 = t371 * qJ(2,3);
t411 = t365 * t371;
t314 = t388 * t349 * t365 + t346 * t408 + (t346 * t411 + (-t358 + 0.1e1) * t349 * t366) * pkin(3);
t355 = 0.1e1 / t365;
t441 = t314 * t355;
t363 = legFrame(2,2);
t347 = sin(t363);
t350 = cos(t363);
t359 = t373 ^ 2;
t368 = sin(qJ(1,2));
t374 = cos(qJ(1,2));
t389 = qJ(2,2) * t368 + t353 * t374;
t407 = t373 * qJ(2,2);
t410 = t367 * t373;
t315 = t389 * t350 * t367 + t347 * t407 + (t347 * t410 + (-t359 + 0.1e1) * t350 * t368) * pkin(3);
t356 = 0.1e1 / t367;
t440 = t315 * t356;
t364 = legFrame(1,2);
t348 = sin(t364);
t351 = cos(t364);
t360 = t375 ^ 2;
t370 = sin(qJ(1,1));
t376 = cos(qJ(1,1));
t390 = qJ(2,1) * t370 + t353 * t376;
t406 = t375 * qJ(2,1);
t409 = t369 * t375;
t316 = t390 * t351 * t369 + t348 * t406 + (t348 * t409 + (-t360 + 0.1e1) * t351 * t370) * pkin(3);
t357 = 0.1e1 / t369;
t439 = t316 * t357;
t317 = (t349 * pkin(3) * t371 - t346 * t388) * t365 + t366 * pkin(3) * (t371 - 0.1e1) * (t371 + 0.1e1) * t346 + t349 * t408;
t438 = t317 * t355;
t318 = (t350 * pkin(3) * t373 - t347 * t389) * t367 + t368 * pkin(3) * (t373 - 0.1e1) * (t373 + 0.1e1) * t347 + t350 * t407;
t437 = t318 * t356;
t319 = (t351 * pkin(3) * t375 - t348 * t390) * t369 + t370 * pkin(3) * (t375 - 0.1e1) * (t375 + 0.1e1) * t348 + t351 * t406;
t436 = t319 * t357;
t340 = t365 * pkin(3) + qJ(2,3);
t326 = t340 * t372 - t353 * t366;
t329 = t442 - (pkin(1) - rSges(2,2)) * m(2);
t337 = 0.1e1 / t340;
t378 = m(2) + m(3);
t320 = (t326 * t378 - t329 * t366) * t337;
t435 = t320 * t355;
t341 = t367 * pkin(3) + qJ(2,2);
t327 = t341 * t374 - t353 * t368;
t338 = 0.1e1 / t341;
t321 = (t327 * t378 - t329 * t368) * t338;
t434 = t321 * t356;
t342 = t369 * pkin(3) + qJ(2,1);
t328 = t342 * t376 - t353 * t370;
t339 = 0.1e1 / t342;
t322 = (t328 * t378 - t329 * t370) * t339;
t433 = t322 * t357;
t432 = t329 * t355;
t431 = t329 * t356;
t430 = t329 * t357;
t429 = t346 * t355;
t428 = t346 * t372;
t427 = t347 * t356;
t426 = t347 * t374;
t425 = t348 * t357;
t424 = t348 * t376;
t423 = t349 * t355;
t422 = t349 * t372;
t421 = t350 * t356;
t420 = t350 * t374;
t419 = t351 * t357;
t418 = t351 * t376;
t417 = t355 * t378;
t386 = 0.1e1 / pkin(3);
t416 = t355 * t386;
t415 = t356 * t378;
t414 = t356 * t386;
t413 = t357 * t378;
t412 = t357 * t386;
t405 = t355 * t445;
t404 = t356 * t444;
t403 = t357 * t443;
t402 = t346 * t416;
t401 = t347 * t414;
t400 = t348 * t412;
t399 = t349 * t416;
t398 = t350 * t414;
t397 = t351 * t412;
t396 = t386 * t405;
t395 = t386 * t404;
t394 = t386 * t403;
t393 = (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + Icges(2,1) + Icges(3,2) + Icges(1,3);
t385 = rSges(3,1) ^ 2;
t387 = pkin(1) ^ 2;
t392 = t385 + t387 + (2 * pkin(1) + t377) * t377;
t391 = rSges(2,3) ^ 2 + t387 + (-2 * pkin(1) + rSges(2,2)) * rSges(2,2);
t383 = rSges(3,2) ^ 2;
t381 = qJ(2,1) ^ 2;
t380 = qJ(2,2) ^ 2;
t379 = qJ(2,3) ^ 2;
t336 = (t383 + t385) * m(3) + Icges(3,3);
t335 = rSges(3,1) * t442 + Icges(3,5);
t334 = -rSges(3,2) * t442 - Icges(3,6);
t333 = (t383 - t385) * m(3) + Icges(3,1) - Icges(3,2);
t325 = t369 * t334 + t335 * t375;
t324 = t367 * t334 + t335 * t373;
t323 = t365 * t334 + t335 * t371;
t313 = (-t325 * t370 + t328 * t443) * t339;
t312 = (-t324 * t368 + t327 * t444) * t338;
t311 = (-t323 * t366 + t326 * t445) * t337;
t310 = t333 * t360 + t409 * t446 + (qJ(2,1) * t447 + t381 + t391) * m(2) + (t381 + 0.2e1 * (rSges(3,1) * t369 + rSges(3,2) * t375) * qJ(2,1) + t392) * m(3) + t393;
t309 = t333 * t359 + t410 * t446 + (qJ(2,2) * t447 + t380 + t391) * m(2) + (t380 + 0.2e1 * (rSges(3,1) * t367 + rSges(3,2) * t373) * qJ(2,2) + t392) * m(3) + t393;
t308 = t333 * t358 + t411 * t446 + (qJ(2,3) * t447 + t379 + t391) * m(2) + (t379 + 0.2e1 * (rSges(3,1) * t365 + rSges(3,2) * t371) * qJ(2,3) + t392) * m(3) + t393;
t307 = -t351 * t394 + (t319 * t413 - t329 * t424) * t339;
t306 = -t350 * t395 + (t318 * t415 - t329 * t426) * t338;
t305 = -t349 * t396 + (t317 * t417 - t329 * t428) * t337;
t304 = -t348 * t394 + (t316 * t413 + t329 * t418) * t339;
t303 = -t347 * t395 + (t315 * t415 + t329 * t420) * t338;
t302 = -t346 * t396 + (t314 * t417 + t329 * t422) * t337;
t301 = (-t310 * t370 + t328 * t329) * t339;
t300 = (-t309 * t368 + t327 * t329) * t338;
t299 = (-t308 * t366 + t326 * t329) * t337;
t298 = -t336 * t397 + (t319 * t403 - t325 * t424) * t339;
t297 = -t336 * t398 + (t318 * t404 - t324 * t426) * t338;
t296 = -t336 * t399 + (t317 * t405 - t323 * t428) * t337;
t295 = -t336 * t400 + (t316 * t403 + t325 * t418) * t339;
t294 = -t336 * t401 + (t315 * t404 + t324 * t420) * t338;
t293 = -t336 * t402 + (t314 * t405 + t323 * t422) * t337;
t292 = -t325 * t397 + (-t310 * t424 + t319 * t430) * t339;
t291 = -t324 * t398 + (-t309 * t426 + t318 * t431) * t338;
t290 = -t323 * t399 + (-t308 * t428 + t317 * t432) * t337;
t289 = -t325 * t400 + (t310 * t418 + t316 * t430) * t339;
t288 = -t324 * t401 + (t309 * t420 + t315 * t431) * t338;
t287 = -t323 * t402 + (t308 * t422 + t314 * t432) * t337;
t1 = [m(4) + (t289 * t418 + t304 * t439) * t339 + (t288 * t420 + t303 * t440) * t338 + (t287 * t422 + t302 * t441) * t337 + (-t293 * t429 - t294 * t427 - t295 * t425) * t386, (-t289 * t424 + t304 * t436) * t339 + (-t288 * t426 + t303 * t437) * t338 + (-t287 * t428 + t302 * t438) * t337 + (-t293 * t423 - t294 * t421 - t295 * t419) * t386, (-t289 * t370 + t304 * t328) * t339 + (-t288 * t368 + t303 * t327) * t338 + (-t287 * t366 + t302 * t326) * t337; (t292 * t418 + t307 * t439) * t339 + (t291 * t420 + t306 * t440) * t338 + (t290 * t422 + t305 * t441) * t337 + (-t296 * t429 - t297 * t427 - t298 * t425) * t386, m(4) + (-t292 * t424 + t307 * t436) * t339 + (-t291 * t426 + t306 * t437) * t338 + (-t290 * t428 + t305 * t438) * t337 + (-t296 * t423 - t297 * t421 - t298 * t419) * t386, (-t292 * t370 + t307 * t328) * t339 + (-t291 * t368 + t306 * t327) * t338 + (-t290 * t366 + t305 * t326) * t337; (t301 * t418 + t316 * t433) * t339 + (t300 * t420 + t315 * t434) * t338 + (t299 * t422 + t314 * t435) * t337 + (-t311 * t429 - t312 * t427 - t313 * t425) * t386, (-t301 * t424 + t319 * t433) * t339 + (-t300 * t426 + t318 * t434) * t338 + (-t299 * t428 + t317 * t435) * t337 + (-t311 * t423 - t312 * t421 - t313 * t419) * t386, m(4) + (-t301 * t370 + t322 * t328) * t339 + (-t300 * t368 + t321 * t327) * t338 + (-t299 * t366 + t320 * t326) * t337;];
MX  = t1;
