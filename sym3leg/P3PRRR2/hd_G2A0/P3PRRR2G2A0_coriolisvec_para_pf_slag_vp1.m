% Calculate vector of centrifugal and coriolis load on the joints for
% P3PRRR2G2P3A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
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
%   pkin=[a3,a4]';
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
% taucX [3x1]
%   forces required to compensate Coriolis and centrifugal joint torques
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:21
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3PRRR2G2P3A0_coriolisvec_para_pf_slag_vp1(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR2G2P3A0_coriolisvec_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRR2G2P3A0_coriolisvec_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR2G2P3A0_coriolisvec_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRR2G2P3A0_coriolisvec_para_pf_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRR2G2P3A0_coriolisvec_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRR2G2P3A0_coriolisvec_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3PRRR2G2P3A0_coriolisvec_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR2G2P3A0_coriolisvec_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR2G2P3A0_coriolisvec_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:21:31
% EndTime: 2020-03-09 21:21:32
% DurationCPUTime: 0.90s
% Computational Cost: add. (2544->126), mult. (4788->239), div. (2610->5), fcn. (5070->21), ass. (0->144)
t402 = legFrame(1,2);
t394 = cos(t402);
t416 = xDP(2);
t464 = t394 * t416;
t391 = sin(t402);
t417 = xDP(1);
t467 = t391 * t417;
t493 = t464 + t467;
t401 = legFrame(2,2);
t393 = cos(t401);
t465 = t393 * t416;
t390 = sin(t401);
t468 = t390 * t417;
t492 = t465 + t468;
t400 = legFrame(3,2);
t392 = cos(t400);
t466 = t392 * t416;
t389 = sin(t400);
t469 = t389 * t417;
t491 = t466 + t469;
t409 = cos(qJ(3,3));
t490 = t409 * pkin(1);
t411 = cos(qJ(3,2));
t489 = t411 * pkin(1);
t413 = cos(qJ(3,1));
t488 = t413 * pkin(1);
t421 = 0.1e1 / pkin(2);
t415 = xDP(3);
t423 = 0.1e1 / pkin(1);
t457 = t415 * t423;
t440 = t421 * t457;
t386 = sin(qJ(2,3) + qJ(3,3));
t404 = sin(qJ(2,3));
t378 = t404 * pkin(1) + pkin(2) * t386;
t403 = sin(qJ(3,3));
t397 = 0.1e1 / t403;
t481 = t378 * t397;
t369 = t440 * t481;
t410 = cos(qJ(2,3));
t460 = t404 * t403;
t372 = (pkin(2) * t409 + pkin(1)) * t410 - pkin(2) * t460;
t484 = t372 * t421;
t426 = t491 * t484;
t375 = t410 * t409 - t460;
t463 = t397 * t423;
t455 = t491 * t375 * t463;
t474 = t386 * t415;
t348 = t369 + (-t426 - t474) * t463 + t455;
t351 = -t426 * t463 + t369;
t487 = t348 * t351;
t387 = sin(qJ(2,2) + qJ(3,2));
t406 = sin(qJ(2,2));
t379 = t406 * pkin(1) + pkin(2) * t387;
t405 = sin(qJ(3,2));
t398 = 0.1e1 / t405;
t480 = t379 * t398;
t370 = t440 * t480;
t412 = cos(qJ(2,2));
t459 = t406 * t405;
t373 = (pkin(2) * t411 + pkin(1)) * t412 - pkin(2) * t459;
t483 = t373 * t421;
t425 = t492 * t483;
t376 = t412 * t411 - t459;
t462 = t398 * t423;
t454 = t492 * t376 * t462;
t472 = t387 * t415;
t349 = t370 + (-t425 - t472) * t462 + t454;
t352 = -t425 * t462 + t370;
t486 = t349 * t352;
t388 = sin(qJ(2,1) + qJ(3,1));
t408 = sin(qJ(2,1));
t380 = t408 * pkin(1) + pkin(2) * t388;
t407 = sin(qJ(3,1));
t399 = 0.1e1 / t407;
t479 = t380 * t399;
t371 = t440 * t479;
t414 = cos(qJ(2,1));
t458 = t408 * t407;
t374 = (pkin(2) * t413 + pkin(1)) * t414 - pkin(2) * t458;
t482 = t374 * t421;
t424 = t493 * t482;
t377 = t414 * t413 - t458;
t461 = t399 * t423;
t453 = t493 * t377 * t461;
t470 = t388 * t415;
t350 = t371 + (-t424 - t470) * t461 + t453;
t353 = -t424 * t461 + t371;
t485 = t350 * t353;
t478 = (t403 * rSges(3,1) + t409 * rSges(3,2)) * t397;
t477 = (t405 * rSges(3,1) + t411 * rSges(3,2)) * t398;
t476 = (t407 * rSges(3,1) + t413 * rSges(3,2)) * t399;
t475 = t386 * t397;
t473 = t387 * t398;
t471 = t388 * t399;
t456 = 0.2e1 * pkin(1) * pkin(2);
t395 = rSges(3,1) ^ 2 + rSges(3,2) ^ 2;
t345 = t369 / 0.2e1 + (-t474 + (-t469 / 0.2e1 - t466 / 0.2e1) * t484) * t463 + t455;
t357 = -t457 * t475 + t455;
t420 = pkin(2) ^ 2;
t422 = pkin(1) ^ 2;
t339 = ((t345 * t409 * t456 + t348 * t420 + t422 * t357) * t421 * t357 + (pkin(2) + t490) * t487) * t463;
t342 = ((-pkin(2) * t348 - t357 * t490) * t357 - pkin(2) * t487) * t463;
t429 = (-rSges(3,1) * t409 + rSges(3,2) * t403) * pkin(1);
t366 = -Icges(3,3) + (-t395 + t429) * m(3);
t385 = t422 + t395;
t433 = -(rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) - Icges(2,3) - Icges(3,3);
t333 = t366 * t339 + (t433 + (-t385 + 0.2e1 * t429) * m(3)) * t342;
t452 = t333 * t375 * t397;
t346 = t370 / 0.2e1 + (-t472 + (-t468 / 0.2e1 - t465 / 0.2e1) * t483) * t462 + t454;
t358 = -t457 * t473 + t454;
t340 = ((t346 * t411 * t456 + t349 * t420 + t422 * t358) * t421 * t358 + (pkin(2) + t489) * t486) * t462;
t343 = ((-pkin(2) * t349 - t358 * t489) * t358 - pkin(2) * t486) * t462;
t428 = (-rSges(3,1) * t411 + rSges(3,2) * t405) * pkin(1);
t367 = -Icges(3,3) + (-t395 + t428) * m(3);
t334 = t367 * t340 + (t433 + (-t385 + 0.2e1 * t428) * m(3)) * t343;
t451 = t334 * t376 * t398;
t347 = t371 / 0.2e1 + (-t470 + (-t467 / 0.2e1 - t464 / 0.2e1) * t482) * t461 + t453;
t359 = -t457 * t471 + t453;
t341 = ((t347 * t413 * t456 + t350 * t420 + t422 * t359) * t421 * t359 + (pkin(2) + t488) * t485) * t461;
t344 = ((-pkin(2) * t350 - t359 * t488) * t359 - pkin(2) * t485) * t461;
t427 = (-rSges(3,1) * t413 + rSges(3,2) * t407) * pkin(1);
t368 = -Icges(3,3) + (-t395 + t427) * m(3);
t335 = t368 * t341 + (t433 + (-t385 + 0.2e1 * t427) * m(3)) * t344;
t450 = t335 * t377 * t399;
t384 = -t395 * m(3) - Icges(3,3);
t336 = t384 * t339 + t366 * t342;
t449 = t336 * t372 * t397;
t337 = t384 * t340 + t367 * t343;
t448 = t337 * t373 * t398;
t338 = t384 * t341 + t368 * t344;
t447 = t338 * t374 * t399;
t446 = t357 ^ 2 * t478;
t445 = t358 ^ 2 * t477;
t444 = t359 ^ 2 * t476;
t439 = t345 * t351 * t478;
t438 = t346 * t352 * t477;
t437 = t347 * t353 * t476;
t436 = t372 * t446;
t435 = t373 * t445;
t434 = t374 * t444;
t432 = -0.2e1 * t375 * t439;
t431 = -0.2e1 * t376 * t438;
t430 = -0.2e1 * t377 * t437;
t1 = [(t389 * t452 + t390 * t451 + t391 * t450 + (-t389 * t449 - t390 * t448 - t391 * t447) * t421) * t423 + (t389 * t432 + t390 * t431 + t391 * t430 + (-t389 * t436 - t390 * t435 - t391 * t434) * t421) * m(3); (t392 * t452 + t393 * t451 + t394 * t450 + (-t392 * t449 - t393 * t448 - t394 * t447) * t421) * t423 + (t392 * t432 + t393 * t431 + t394 * t430 + (-t392 * t436 - t393 * t435 - t394 * t434) * t421) * m(3); (-t333 * t475 - t334 * t473 - t335 * t471 + (t336 * t481 + t337 * t480 + t338 * t479) * t421) * t423 + (0.2e1 * t386 * t439 + 0.2e1 * t387 * t438 + 0.2e1 * t388 * t437 + (t378 * t446 + t379 * t445 + t380 * t444) * t421) * m(3);];
taucX  = t1;
