% Calculate inertia matrix for parallel robot
% P3RRPRR12V1G3A0
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
% Datum: 2020-08-06 19:11
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RRPRR12V1G3A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(4,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V1G3A0_inertia_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V1G3A0_inertia_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RRPRR12V1G3A0_inertia_para_pf_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR12V1G3A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRPRR12V1G3A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RRPRR12V1G3A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V1G3A0_inertia_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V1G3A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:09:47
% EndTime: 2020-08-06 19:09:49
% DurationCPUTime: 1.95s
% Computational Cost: add. (3387->244), mult. (4428->407), div. (432->6), fcn. (2655->18), ass. (0->183)
t480 = 0.1e1 / qJ(3,1);
t465 = sin(qJ(2,1));
t471 = cos(qJ(2,1));
t474 = pkin(1) + pkin(2);
t434 = t465 * qJ(3,1) + t474 * t471;
t466 = sin(qJ(1,1));
t472 = cos(qJ(1,1));
t571 = t472 * pkin(4) + t434 * t466;
t533 = t571 * t480;
t478 = 0.1e1 / qJ(3,2);
t463 = sin(qJ(2,2));
t469 = cos(qJ(2,2));
t433 = t463 * qJ(3,2) + t474 * t469;
t464 = sin(qJ(1,2));
t470 = cos(qJ(1,2));
t573 = t470 * pkin(4) + t433 * t464;
t534 = t573 * t478;
t476 = 0.1e1 / qJ(3,3);
t461 = sin(qJ(2,3));
t467 = cos(qJ(2,3));
t432 = t461 * qJ(3,3) + t474 * t467;
t462 = sin(qJ(1,3));
t468 = cos(qJ(1,3));
t575 = t468 * pkin(4) + t432 * t462;
t535 = t575 * t476;
t562 = t462 * pkin(4);
t574 = t432 * t468 - t562;
t561 = t464 * pkin(4);
t572 = t433 * t470 - t561;
t560 = t466 * pkin(4);
t570 = t434 * t472 - t560;
t566 = 2 * rSges(3,3);
t569 = (t566 + qJ(3,3)) * qJ(3,3);
t568 = (t566 + qJ(3,2)) * qJ(3,2);
t567 = (t566 + qJ(3,1)) * qJ(3,1);
t565 = m(2) * rSges(2,3);
t564 = rSges(3,2) * m(3);
t473 = pkin(1) + rSges(3,1);
t563 = m(3) * t473;
t556 = rSges(3,2) * t461;
t555 = rSges(3,2) * t463;
t554 = rSges(3,2) * t465;
t520 = t461 * t468;
t423 = qJ(3,3) * t520 - t562;
t458 = legFrame(3,2);
t442 = sin(t458);
t445 = cos(t458);
t452 = t467 ^ 2;
t511 = t474 * t461;
t517 = t468 * t474;
t532 = t442 * qJ(3,3);
t398 = (t445 * t517 - t532) * t452 + (t423 * t445 + t442 * t511) * t467 + t532;
t553 = t398 * t476;
t519 = t463 * t470;
t424 = qJ(3,2) * t519 - t561;
t459 = legFrame(2,2);
t443 = sin(t459);
t446 = cos(t459);
t453 = t469 ^ 2;
t510 = t474 * t463;
t516 = t470 * t474;
t530 = t443 * qJ(3,2);
t399 = (t446 * t516 - t530) * t453 + (t424 * t446 + t443 * t510) * t469 + t530;
t552 = t399 * t478;
t518 = t465 * t472;
t425 = qJ(3,1) * t518 - t560;
t460 = legFrame(1,2);
t444 = sin(t460);
t447 = cos(t460);
t454 = t471 ^ 2;
t509 = t474 * t465;
t515 = t472 * t474;
t528 = t444 * qJ(3,1);
t400 = (t447 * t515 - t528) * t454 + (t425 * t447 + t444 * t509) * t471 + t528;
t551 = t400 * t480;
t526 = t445 * qJ(3,3);
t401 = (-t442 * t517 - t526) * t452 + (-t442 * t423 + t445 * t511) * t467 + t526;
t550 = t401 * t476;
t524 = t446 * qJ(3,2);
t402 = (-t443 * t516 - t524) * t453 + (-t443 * t424 + t446 * t510) * t469 + t524;
t549 = t402 * t478;
t522 = t447 * qJ(3,1);
t403 = (-t444 * t515 - t522) * t454 + (-t425 * t444 + t447 * t509) * t471 + t522;
t548 = t403 * t480;
t429 = -qJ(3,3) * t467 + t511;
t404 = t429 * t442 + t445 * t574;
t547 = t404 * t476;
t430 = -qJ(3,2) * t469 + t510;
t405 = t430 * t443 + t446 * t572;
t546 = t405 * t478;
t431 = -qJ(3,1) * t471 + t509;
t406 = t431 * t444 + t447 * t570;
t545 = t406 * t480;
t407 = t429 * t445 - t442 * t574;
t544 = t407 * t476;
t408 = t430 * t446 - t443 * t572;
t543 = t408 * t478;
t409 = t431 * t447 - t444 * t570;
t542 = t409 * t480;
t422 = rSges(2,1) * t565 + rSges(3,2) * t563 - Icges(3,4) - Icges(2,5);
t455 = rSges(3,3) + qJ(3,3);
t487 = -rSges(2,2) * t565 + Icges(2,6) - Icges(3,6);
t410 = (t455 * t564 + t487) * t467 - t422 * t461;
t541 = t410 * t476;
t456 = rSges(3,3) + qJ(3,2);
t411 = (t456 * t564 + t487) * t469 - t422 * t463;
t540 = t411 * t478;
t457 = rSges(3,3) + qJ(3,1);
t412 = (t457 * t564 + t487) * t471 - t422 * t465;
t539 = t412 * t480;
t481 = rSges(3,3) ^ 2;
t491 = pkin(1) ^ 2 + t481 + (2 * pkin(1) + rSges(3,1)) * rSges(3,1);
t483 = rSges(2,2) ^ 2;
t485 = rSges(2,1) ^ 2;
t495 = Icges(3,2) + Icges(2,3) + (t483 + t485) * m(2);
t413 = (t491 + t569) * m(3) + t495;
t538 = t413 * t476;
t414 = (t491 + t568) * m(3) + t495;
t537 = t414 * t478;
t415 = (t491 + t567) * m(3) + t495;
t536 = t415 * t480;
t531 = t442 * t462;
t529 = t443 * t464;
t527 = t444 * t466;
t525 = t445 * t462;
t523 = t446 * t464;
t521 = t447 * t466;
t514 = t473 * t476;
t513 = t473 * t478;
t512 = t473 * t480;
t508 = -Icges(2,1) - Icges(3,1);
t507 = rSges(3,2) ^ 2 + t481;
t506 = rSges(3,3) - t473;
t505 = rSges(3,3) + t473;
t504 = m(3) * t514;
t503 = m(3) * t513;
t502 = m(3) * t512;
t501 = t462 * t556;
t500 = t464 * t555;
t499 = t466 * t554;
t498 = t467 * t535;
t497 = t469 * t534;
t496 = t471 * t533;
t494 = m(3) * t476 * t556;
t493 = m(3) * t478 * t555;
t492 = m(3) * t480 * t554;
t490 = Icges(1,3) + (rSges(2,3) ^ 2 + t483) * m(2) + (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) - t508;
t489 = (-t483 + t485) * m(2) + Icges(2,2) + Icges(3,3) + t508;
t488 = -m(2) * rSges(2,1) * rSges(2,2) + Icges(2,4) - Icges(3,5);
t428 = 0.1e1 / t434;
t427 = 0.1e1 / t433;
t426 = 0.1e1 / t432;
t397 = (-(qJ(3,1) + t505) * (qJ(3,1) + t506) * m(3) + t489) * t454 + 0.2e1 * (t457 * t563 + t488) * t465 * t471 + (t507 + t567) * m(3) + t490;
t396 = (-(qJ(3,2) + t505) * (qJ(3,2) + t506) * m(3) + t489) * t453 + 0.2e1 * (t456 * t563 + t488) * t463 * t469 + (t507 + t568) * m(3) + t490;
t395 = (-(qJ(3,3) + t505) * (qJ(3,3) + t506) * m(3) + t489) * t452 + 0.2e1 * (t455 * t563 + t488) * t461 * t467 + (t507 + t569) * m(3) + t490;
t394 = (-t533 + (-rSges(3,2) * t518 + t473 * t496) * t428) * m(3);
t393 = (-t534 + (-rSges(3,2) * t519 + t473 * t497) * t427) * m(3);
t392 = (-t535 + (-rSges(3,2) * t520 + t473 * t498) * t426) * m(3);
t391 = t571 * t502 + (-t412 * t472 - t415 * t496) * t428;
t390 = t573 * t503 + (-t411 * t470 - t414 * t497) * t427;
t389 = t575 * t504 + (-t410 * t468 - t413 * t498) * t426;
t388 = (t545 + (-t400 * t512 - t447 * t499) * t428) * m(3);
t387 = (t546 + (-t399 * t513 - t446 * t500) * t427) * m(3);
t386 = (t547 + (-t398 * t514 - t445 * t501) * t426) * m(3);
t385 = (t542 + (-t403 * t512 + t444 * t499) * t428) * m(3);
t384 = (t543 + (-t402 * t513 + t443 * t500) * t427) * m(3);
t383 = (t544 + (-t401 * t514 + t442 * t501) * t426) * m(3);
t382 = -t571 * t492 + (-t397 * t472 - t412 * t496) * t428;
t381 = -t573 * t493 + (-t396 * t470 - t411 * t497) * t427;
t380 = -t575 * t494 + (-t395 * t468 - t410 * t498) * t426;
t379 = -t406 * t502 + (t400 * t536 - t412 * t521) * t428;
t378 = -t405 * t503 + (t399 * t537 - t411 * t523) * t427;
t377 = -t404 * t504 + (t398 * t538 - t410 * t525) * t426;
t376 = -t409 * t502 + (t403 * t536 + t412 * t527) * t428;
t375 = -t408 * t503 + (t402 * t537 + t411 * t529) * t427;
t374 = -t407 * t504 + (t401 * t538 + t410 * t531) * t426;
t373 = t406 * t492 + (-t397 * t521 + t400 * t539) * t428;
t372 = t405 * t493 + (-t396 * t523 + t399 * t540) * t427;
t371 = t404 * t494 + (-t395 * t525 + t398 * t541) * t426;
t370 = t409 * t492 + (t397 * t527 + t403 * t539) * t428;
t369 = t408 * t493 + (t396 * t529 + t402 * t540) * t427;
t368 = t407 * t494 + (t395 * t531 + t401 * t541) * t426;
t1 = [t386 * t547 + t387 * t546 + t388 * t545 + m(4) + (-t373 * t521 + t379 * t551) * t428 + (-t372 * t523 + t378 * t552) * t427 + (-t371 * t525 + t377 * t553) * t426, t386 * t544 + t387 * t543 + t388 * t542 + (t373 * t527 + t379 * t548) * t428 + (t372 * t529 + t378 * t549) * t427 + (t371 * t531 + t377 * t550) * t426, -t386 * t535 - t387 * t534 - t388 * t533 + (-t373 * t472 - t379 * t496) * t428 + (-t372 * t470 - t378 * t497) * t427 + (-t371 * t468 - t377 * t498) * t426; t383 * t547 + t384 * t546 + t385 * t545 + (-t370 * t521 + t376 * t551) * t428 + (-t369 * t523 + t375 * t552) * t427 + (-t368 * t525 + t374 * t553) * t426, t383 * t544 + t384 * t543 + t385 * t542 + m(4) + (t370 * t527 + t376 * t548) * t428 + (t369 * t529 + t375 * t549) * t427 + (t368 * t531 + t374 * t550) * t426, -t383 * t535 - t384 * t534 - t385 * t533 + (-t370 * t472 - t376 * t496) * t428 + (-t369 * t470 - t375 * t497) * t427 + (-t368 * t468 - t374 * t498) * t426; t392 * t547 + t393 * t546 + t394 * t545 + (-t382 * t521 + t391 * t551) * t428 + (-t381 * t523 + t390 * t552) * t427 + (-t380 * t525 + t389 * t553) * t426, t392 * t544 + t393 * t543 + t394 * t542 + (t382 * t527 + t391 * t548) * t428 + (t381 * t529 + t390 * t549) * t427 + (t380 * t531 + t389 * t550) * t426, -t392 * t535 - t393 * t534 - t394 * t533 + m(4) + (-t382 * t472 - t391 * t496) * t428 + (-t381 * t470 - t390 * t497) * t427 + (-t380 * t468 - t389 * t498) * t426;];
MX  = t1;
