% Calculate inertia matrix for parallel robot
% P3RPRRR9V1G1A0
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
% Datum: 2020-08-06 18:48
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RPRRR9V1G1A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR9V1G1A0_inertia_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR9V1G1A0_inertia_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR9V1G1A0_inertia_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR9V1G1A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RPRRR9V1G1A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RPRRR9V1G1A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR9V1G1A0_inertia_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR9V1G1A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:46:49
% EndTime: 2020-08-06 18:46:50
% DurationCPUTime: 0.98s
% Computational Cost: add. (3339->209), mult. (4323->329), div. (270->7), fcn. (1734->40), ass. (0->147)
t568 = 2 * pkin(1);
t567 = Icges(3,2) / 0.2e1;
t497 = (pkin(7) + qJ(3,3));
t461 = 2 * t497;
t456 = sin(t461);
t473 = sin(t497);
t521 = 2 * pkin(7);
t494 = t521 + qJ(3,3);
t538 = sin(t494) + sin(qJ(3,3));
t434 = t538 * pkin(2) + pkin(3) * t456 + t473 * t568;
t562 = t434 / 0.2e1;
t498 = pkin(7) + qJ(3,2);
t462 = 2 * t498;
t457 = sin(t462);
t474 = sin(t498);
t495 = t521 + qJ(3,2);
t537 = sin(t495) + sin(qJ(3,2));
t435 = t537 * pkin(2) + pkin(3) * t457 + t474 * t568;
t561 = t435 / 0.2e1;
t499 = pkin(7) + qJ(3,1);
t463 = 2 * t499;
t458 = sin(t463);
t475 = sin(t499);
t496 = t521 + qJ(3,1);
t536 = sin(t496) + sin(qJ(3,1));
t436 = t536 * pkin(2) + pkin(3) * t458 + t475 * t568;
t560 = t436 / 0.2e1;
t565 = 4 * rSges(2,3);
t564 = m(2) / 0.2e1;
t563 = rSges(2,2) * m(2);
t522 = rSges(3,2) ^ 2;
t524 = rSges(3,1) ^ 2;
t559 = (-t522 + t524) * m(3) / 0.2e1 - Icges(3,1) / 0.2e1 + t567;
t518 = m(2) + m(3);
t558 = t518 / 0.2e1;
t549 = pkin(5) + qJ(2,3);
t488 = rSges(3,3) + t549;
t557 = m(3) * t488;
t550 = pkin(5) + qJ(2,2);
t489 = rSges(3,3) + t550;
t556 = m(3) * t489;
t551 = pkin(5) + qJ(2,1);
t490 = rSges(3,3) + t551;
t555 = m(3) * t490;
t476 = cos(t497);
t554 = pkin(3) * t476;
t477 = cos(t498);
t553 = pkin(3) * t477;
t478 = cos(t499);
t552 = pkin(3) * t478;
t526 = 0.1e1 / pkin(3);
t548 = ((t522 + t524) * m(3) + Icges(3,3)) * t526;
t501 = cos(pkin(7));
t452 = (m(2) * rSges(2,1) + pkin(2) * m(3)) * t501;
t467 = 0.1e1 / t476;
t491 = -pkin(6) - t549;
t485 = 0.1e1 / t491;
t547 = t467 * t485;
t468 = 0.1e1 / t477;
t492 = -pkin(6) - t550;
t486 = 0.1e1 / t492;
t546 = t468 * t486;
t469 = 0.1e1 / t478;
t493 = -pkin(6) - t551;
t487 = 0.1e1 / t493;
t545 = t469 * t487;
t544 = t473 * t485;
t543 = t474 * t486;
t542 = t475 * t487;
t431 = (-rSges(3,2) * t557 + Icges(3,6)) * t476 - (rSges(3,1) * t557 - Icges(3,5)) * t473;
t541 = t526 * t431;
t432 = (-rSges(3,2) * t556 + Icges(3,6)) * t477 - (rSges(3,1) * t556 - Icges(3,5)) * t474;
t540 = t526 * t432;
t433 = (-rSges(3,2) * t555 + Icges(3,6)) * t478 - (rSges(3,1) * t555 - Icges(3,5)) * t475;
t539 = t526 * t433;
t466 = sin(pkin(7)) * t563;
t520 = 2 * pkin(1) ^ 2;
t523 = rSges(2,2) ^ 2;
t525 = rSges(2,1) ^ 2;
t534 = (2 * rSges(2,3) ^ 2) + t520 + t523 + t525;
t533 = -pkin(1) * t518 - t452 + t466;
t532 = rSges(3,1) * t476 - rSges(3,2) * t473;
t531 = rSges(3,1) * t477 - rSges(3,2) * t474;
t530 = rSges(3,1) * t478 - rSges(3,2) * t475;
t527 = pkin(2) ^ 2;
t529 = t520 / 0.2e1 + t522 / 0.2e1 + t524 / 0.2e1 + t527 / 0.2e1;
t528 = Icges(1,3) + (m(3) * t527 + (-t523 + t525) * m(2) - Icges(2,1) + Icges(2,2)) * cos(t521) / 0.2e1 + (-rSges(2,1) * t563 + Icges(2,4)) * sin(t521) + t452 * t568 + (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) - 0.2e1 * pkin(1) * t466 + t567 + Icges(2,2) / 0.2e1 + Icges(3,1) / 0.2e1 + Icges(2,1) / 0.2e1;
t513 = cos(qJ(1,1));
t512 = cos(qJ(1,2));
t511 = cos(qJ(1,3));
t510 = sin(qJ(1,1));
t508 = sin(qJ(1,2));
t506 = sin(qJ(1,3));
t504 = legFrame(1,3);
t503 = legFrame(2,3);
t502 = legFrame(3,3);
t484 = cos(t504);
t483 = cos(t503);
t482 = cos(t502);
t481 = sin(t504);
t480 = sin(t503);
t479 = sin(t502);
t465 = -m(3) * rSges(3,1) * rSges(3,2) + Icges(3,4);
t460 = t501 * pkin(2) + pkin(1);
t449 = t481 * t513 + t484 * t510;
t448 = t480 * t512 + t483 * t508;
t447 = t479 * t511 + t482 * t506;
t446 = -t481 * t510 + t484 * t513;
t445 = -t480 * t508 + t483 * t512;
t444 = -t479 * t506 + t482 * t511;
t442 = t460 * t513 - t510 * t493;
t441 = t460 * t512 - t508 * t492;
t440 = t460 * t511 - t506 * t491;
t439 = t510 * t460 + t513 * t493;
t438 = t508 * t460 + t512 * t492;
t437 = t506 * t460 + t511 * t491;
t430 = -t530 * m(3) + t533;
t429 = -t531 * m(3) + t533;
t428 = -t532 * m(3) + t533;
t427 = t439 * t484 + t442 * t481 + t449 * t552;
t426 = t438 * t483 + t441 * t480 + t448 * t553;
t425 = t437 * t482 + t440 * t479 + t447 * t554;
t424 = -t439 * t481 + t442 * t484 + t446 * t552;
t423 = -t438 * t480 + t441 * t483 + t445 * t553;
t422 = -t437 * t479 + t440 * t482 + t444 * t554;
t421 = (t475 * t430 + t436 * t558) * t545;
t420 = (t474 * t429 + t435 * t558) * t546;
t419 = (t473 * t428 + t434 * t558) * t547;
t418 = (t427 * t518 + t430 * t449) * t487;
t417 = (t426 * t518 + t429 * t448) * t486;
t416 = (t425 * t518 + t428 * t447) * t485;
t415 = (t424 * t518 + t430 * t446) * t487;
t414 = (t423 * t518 + t429 * t445) * t486;
t413 = (t422 * t518 + t428 * t444) * t485;
t412 = t528 + t465 * t458 + (t490 ^ 2 + t530 * t568 + (-t536 * rSges(3,2) + (cos(t496) + cos(qJ(3,1))) * rSges(3,1)) * pkin(2) + t529) * m(3) + (((t565 + 2 * qJ(2,1)) * qJ(2,1)) + t534) * t564 + cos(t463) * t559;
t411 = t528 + t465 * t457 + (t489 ^ 2 + t531 * t568 + (-t537 * rSges(3,2) + (cos(qJ(3,2)) + cos(t495)) * rSges(3,1)) * pkin(2) + t529) * m(3) + (((t565 + 2 * qJ(2,2)) * qJ(2,2)) + t534) * t564 + cos(t462) * t559;
t410 = t528 + t465 * t456 + (t488 ^ 2 + t532 * t568 + (-t538 * rSges(3,2) + (cos(t494) + cos(qJ(3,3))) * rSges(3,1)) * pkin(2) + t529) * m(3) + (((t565 + 2 * qJ(2,3)) * qJ(2,3)) + t534) * t564 + cos(t461) * t559;
t409 = (t539 - (t475 * t412 + t430 * t560) * t487) * t469;
t408 = (t540 - (t474 * t411 + t429 * t561) * t486) * t468;
t407 = (t541 - (t473 * t410 + t428 * t562) * t485) * t467;
t406 = (t412 * t449 + t427 * t430) * t487;
t405 = (t411 * t448 + t426 * t429) * t486;
t404 = (t410 * t447 + t425 * t428) * t485;
t403 = (t412 * t446 + t424 * t430) * t487;
t402 = (t411 * t445 + t423 * t429) * t486;
t401 = (t410 * t444 + t422 * t428) * t485;
t1 = [m(4) - (-t403 * t446 - t415 * t424) * t487 - (-t402 * t445 - t414 * t423) * t486 - (-t401 * t444 - t413 * t422) * t485, -(-t403 * t449 - t415 * t427) * t487 - (-t402 * t448 - t414 * t426) * t486 - (-t401 * t447 - t413 * t425) * t485, -(-t403 * t475 - t415 * t560 + t446 * t539) * t545 - (-t402 * t474 - t414 * t561 + t445 * t540) * t546 - (-t401 * t473 - t413 * t562 + t444 * t541) * t547; -(-t406 * t446 - t418 * t424) * t487 - (-t405 * t445 - t417 * t423) * t486 - (-t404 * t444 - t416 * t422) * t485, m(4) - (-t406 * t449 - t418 * t427) * t487 - (-t405 * t448 - t417 * t426) * t486 - (-t404 * t447 - t416 * t425) * t485, -(-t406 * t475 - t418 * t560 + t449 * t539) * t545 - (-t405 * t474 - t417 * t561 + t448 * t540) * t546 - (-t404 * t473 - t416 * t562 + t447 * t541) * t547; -(t409 * t446 - t421 * t424) * t487 - (t408 * t445 - t420 * t423) * t486 - (t407 * t444 - t419 * t422) * t485, -(t409 * t449 - t421 * t427) * t487 - (t408 * t448 - t420 * t426) * t486 - (t407 * t447 - t419 * t425) * t485, m(4) + (-t409 * t542 + t421 * t487 * t560 + (-t433 * t542 + t548) * t526 * t469) * t469 + (-t408 * t543 + t420 * t486 * t561 + (-t432 * t543 + t548) * t526 * t468) * t468 + (-t407 * t544 + t419 * t485 * t562 + (-t431 * t544 + t548) * t526 * t467) * t467;];
MX  = t1;
