% Calculate inertia matrix for parallel robot
% P3RRPRR12V1G2A0
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d1,d4]';
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
% Datum: 2020-08-06 19:07
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RRPRR12V1G2A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(4,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V1G2A0_inertia_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V1G2A0_inertia_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RRPRR12V1G2A0_inertia_para_pf_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR12V1G2A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRPRR12V1G2A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RRPRR12V1G2A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V1G2A0_inertia_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V1G2A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:05:31
% EndTime: 2020-08-06 19:05:32
% DurationCPUTime: 1.45s
% Computational Cost: add. (3387->244), mult. (4500->407), div. (432->6), fcn. (2727->18), ass. (0->183)
t478 = 0.1e1 / qJ(3,1);
t463 = sin(qJ(2,1));
t469 = cos(qJ(2,1));
t472 = pkin(1) + pkin(2);
t435 = t463 * qJ(3,1) + t472 * t469;
t464 = sin(qJ(1,1));
t470 = cos(qJ(1,1));
t489 = -t464 * pkin(4) + t435 * t470;
t534 = t489 * t478;
t476 = 0.1e1 / qJ(3,2);
t461 = sin(qJ(2,2));
t467 = cos(qJ(2,2));
t434 = t461 * qJ(3,2) + t472 * t467;
t462 = sin(qJ(1,2));
t468 = cos(qJ(1,2));
t490 = -t462 * pkin(4) + t434 * t468;
t535 = t490 * t476;
t474 = 0.1e1 / qJ(3,3);
t459 = sin(qJ(2,3));
t465 = cos(qJ(2,3));
t433 = t459 * qJ(3,3) + t472 * t465;
t460 = sin(qJ(1,3));
t466 = cos(qJ(1,3));
t491 = -t460 * pkin(4) + t433 * t466;
t536 = t491 * t474;
t564 = 2 * rSges(3,3);
t567 = (t564 + qJ(3,3)) * qJ(3,3);
t566 = (t564 + qJ(3,2)) * qJ(3,2);
t565 = (t564 + qJ(3,1)) * qJ(3,1);
t563 = m(2) * rSges(2,3);
t562 = rSges(3,2) * m(3);
t471 = pkin(1) + rSges(3,1);
t561 = m(3) * t471;
t560 = t466 * pkin(4);
t559 = t468 * pkin(4);
t558 = t470 * pkin(4);
t557 = rSges(3,2) * t459;
t556 = rSges(3,2) * t461;
t555 = rSges(3,2) * t463;
t521 = t459 * t460;
t424 = qJ(3,3) * t521 + t560;
t456 = legFrame(3,2);
t440 = sin(t456);
t443 = cos(t456);
t450 = t465 ^ 2;
t512 = t472 * t459;
t520 = t460 * t472;
t533 = t440 * qJ(3,3);
t396 = (t443 * t520 - t533) * t450 + (t424 * t443 + t440 * t512) * t465 + t533;
t554 = t396 * t474;
t519 = t461 * t462;
t425 = qJ(3,2) * t519 + t559;
t457 = legFrame(2,2);
t441 = sin(t457);
t444 = cos(t457);
t451 = t467 ^ 2;
t511 = t472 * t461;
t518 = t462 * t472;
t531 = t441 * qJ(3,2);
t397 = (t444 * t518 - t531) * t451 + (t425 * t444 + t441 * t511) * t467 + t531;
t553 = t397 * t476;
t517 = t463 * t464;
t426 = qJ(3,1) * t517 + t558;
t458 = legFrame(1,2);
t442 = sin(t458);
t445 = cos(t458);
t452 = t469 ^ 2;
t510 = t472 * t463;
t516 = t464 * t472;
t529 = t442 * qJ(3,1);
t398 = (t445 * t516 - t529) * t452 + (t426 * t445 + t442 * t510) * t469 + t529;
t552 = t398 * t478;
t527 = t443 * qJ(3,3);
t399 = (-t440 * t520 - t527) * t450 + (-t440 * t424 + t443 * t512) * t465 + t527;
t551 = t399 * t474;
t525 = t444 * qJ(3,2);
t400 = (-t441 * t518 - t525) * t451 + (-t441 * t425 + t444 * t511) * t467 + t525;
t550 = t400 * t476;
t523 = t445 * qJ(3,1);
t401 = (-t442 * t516 - t523) * t452 + (-t442 * t426 + t445 * t510) * t469 + t523;
t549 = t401 * t478;
t414 = t433 * t460 + t560;
t430 = -t465 * qJ(3,3) + t512;
t402 = t414 * t443 + t440 * t430;
t548 = t402 * t474;
t403 = -t414 * t440 + t430 * t443;
t547 = t403 * t474;
t415 = t434 * t462 + t559;
t431 = -t467 * qJ(3,2) + t511;
t404 = t415 * t444 + t441 * t431;
t546 = t404 * t476;
t405 = -t415 * t441 + t431 * t444;
t545 = t405 * t476;
t416 = t435 * t464 + t558;
t432 = -t469 * qJ(3,1) + t510;
t406 = t416 * t445 + t442 * t432;
t544 = t406 * t478;
t407 = -t416 * t442 + t432 * t445;
t543 = t407 * t478;
t423 = rSges(2,1) * t563 + rSges(3,2) * t561 - Icges(3,4) - Icges(2,5);
t453 = rSges(3,3) + qJ(3,3);
t485 = -rSges(2,2) * t563 + Icges(2,6) - Icges(3,6);
t408 = (t453 * t562 + t485) * t465 - t423 * t459;
t542 = t408 * t474;
t454 = rSges(3,3) + qJ(3,2);
t409 = (t454 * t562 + t485) * t467 - t423 * t461;
t541 = t409 * t476;
t455 = rSges(3,3) + qJ(3,1);
t410 = (t455 * t562 + t485) * t469 - t423 * t463;
t540 = t410 * t478;
t479 = rSges(3,3) ^ 2;
t492 = pkin(1) ^ 2 + t479 + (2 * pkin(1) + rSges(3,1)) * rSges(3,1);
t481 = rSges(2,2) ^ 2;
t483 = rSges(2,1) ^ 2;
t496 = Icges(3,2) + Icges(2,3) + (t481 + t483) * m(2);
t411 = (t492 + t567) * m(3) + t496;
t539 = t411 * t474;
t412 = (t492 + t566) * m(3) + t496;
t538 = t412 * t476;
t413 = (t492 + t565) * m(3) + t496;
t537 = t413 * t478;
t532 = t440 * t466;
t530 = t441 * t468;
t528 = t442 * t470;
t526 = t443 * t466;
t524 = t444 * t468;
t522 = t445 * t470;
t515 = t471 * t474;
t514 = t471 * t476;
t513 = t471 * t478;
t509 = -Icges(2,1) - Icges(3,1);
t508 = rSges(3,2) ^ 2 + t479;
t507 = rSges(3,3) - t471;
t506 = rSges(3,3) + t471;
t505 = m(3) * t515;
t504 = m(3) * t514;
t503 = m(3) * t513;
t502 = t466 * t557;
t501 = t468 * t556;
t500 = t470 * t555;
t499 = t465 * t536;
t498 = t467 * t535;
t497 = t469 * t534;
t495 = m(3) * t474 * t557;
t494 = m(3) * t476 * t556;
t493 = m(3) * t478 * t555;
t488 = Icges(1,3) + (rSges(2,3) ^ 2 + t481) * m(2) + (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) - t509;
t487 = (-t481 + t483) * m(2) + Icges(2,2) + Icges(3,3) + t509;
t486 = -m(2) * rSges(2,1) * rSges(2,2) + Icges(2,4) - Icges(3,5);
t429 = 0.1e1 / t435;
t428 = 0.1e1 / t434;
t427 = 0.1e1 / t433;
t395 = (-(qJ(3,1) + t506) * (qJ(3,1) + t507) * m(3) + t487) * t452 + 0.2e1 * t463 * (t455 * t561 + t486) * t469 + (t508 + t565) * m(3) + t488;
t394 = (-(qJ(3,2) + t506) * (qJ(3,2) + t507) * m(3) + t487) * t451 + 0.2e1 * t461 * (t454 * t561 + t486) * t467 + (t508 + t566) * m(3) + t488;
t393 = (-(qJ(3,3) + t506) * (qJ(3,3) + t507) * m(3) + t487) * t450 + 0.2e1 * t459 * (t453 * t561 + t486) * t465 + (t508 + t567) * m(3) + t488;
t392 = (t534 + (-rSges(3,2) * t517 - t471 * t497) * t429) * m(3);
t391 = (t535 + (-rSges(3,2) * t519 - t471 * t498) * t428) * m(3);
t390 = (t536 + (-rSges(3,2) * t521 - t471 * t499) * t427) * m(3);
t389 = -t489 * t503 + (-t410 * t464 + t413 * t497) * t429;
t388 = -t490 * t504 + (-t409 * t462 + t412 * t498) * t428;
t387 = -t491 * t505 + (-t408 * t460 + t411 * t499) * t427;
t386 = (t544 + (-t398 * t513 + t445 * t500) * t429) * m(3);
t385 = (t546 + (-t397 * t514 + t444 * t501) * t428) * m(3);
t384 = (t548 + (-t396 * t515 + t443 * t502) * t427) * m(3);
t383 = (t543 + (-t401 * t513 - t442 * t500) * t429) * m(3);
t382 = (t545 + (-t400 * t514 - t441 * t501) * t428) * m(3);
t381 = (t547 + (-t399 * t515 - t440 * t502) * t427) * m(3);
t380 = t489 * t493 + (-t395 * t464 + t410 * t497) * t429;
t379 = t490 * t494 + (-t394 * t462 + t409 * t498) * t428;
t378 = t491 * t495 + (-t393 * t460 + t408 * t499) * t427;
t377 = -t406 * t503 + (t398 * t537 + t410 * t522) * t429;
t376 = -t404 * t504 + (t397 * t538 + t409 * t524) * t428;
t375 = -t402 * t505 + (t396 * t539 + t408 * t526) * t427;
t374 = -t407 * t503 + (t401 * t537 - t410 * t528) * t429;
t373 = -t405 * t504 + (t400 * t538 - t409 * t530) * t428;
t372 = -t403 * t505 + (t399 * t539 - t408 * t532) * t427;
t371 = t406 * t493 + (t395 * t522 + t398 * t540) * t429;
t370 = t404 * t494 + (t394 * t524 + t397 * t541) * t428;
t369 = t402 * t495 + (t393 * t526 + t396 * t542) * t427;
t368 = t407 * t493 + (-t395 * t528 + t401 * t540) * t429;
t367 = t405 * t494 + (-t394 * t530 + t400 * t541) * t428;
t366 = t403 * t495 + (-t393 * t532 + t399 * t542) * t427;
t1 = [t384 * t548 + t385 * t546 + t386 * t544 + m(4) + (t371 * t522 + t377 * t552) * t429 + (t370 * t524 + t376 * t553) * t428 + (t369 * t526 + t375 * t554) * t427, t384 * t547 + t385 * t545 + t386 * t543 + (-t371 * t528 + t377 * t549) * t429 + (-t370 * t530 + t376 * t550) * t428 + (-t369 * t532 + t375 * t551) * t427, t384 * t536 + t385 * t535 + t386 * t534 + (-t371 * t464 + t377 * t497) * t429 + (-t370 * t462 + t376 * t498) * t428 + (-t369 * t460 + t375 * t499) * t427; t381 * t548 + t382 * t546 + t383 * t544 + (t368 * t522 + t374 * t552) * t429 + (t367 * t524 + t373 * t553) * t428 + (t366 * t526 + t372 * t554) * t427, t381 * t547 + t382 * t545 + t383 * t543 + m(4) + (-t368 * t528 + t374 * t549) * t429 + (-t367 * t530 + t373 * t550) * t428 + (-t366 * t532 + t372 * t551) * t427, t381 * t536 + t382 * t535 + t383 * t534 + (-t368 * t464 + t374 * t497) * t429 + (-t367 * t462 + t373 * t498) * t428 + (-t366 * t460 + t372 * t499) * t427; t390 * t548 + t391 * t546 + t392 * t544 + (t380 * t522 + t389 * t552) * t429 + (t379 * t524 + t388 * t553) * t428 + (t378 * t526 + t387 * t554) * t427, t390 * t547 + t391 * t545 + t392 * t543 + (-t380 * t528 + t389 * t549) * t429 + (-t379 * t530 + t388 * t550) * t428 + (-t378 * t532 + t387 * t551) * t427, t390 * t536 + t391 * t535 + t392 * t534 + m(4) + (-t380 * t464 + t389 * t497) * t429 + (-t379 * t462 + t388 * t498) * t428 + (-t378 * t460 + t387 * t499) * t427;];
MX  = t1;
