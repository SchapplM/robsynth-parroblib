% Calculate inertia matrix for parallel robot
% P3RRPRR12V1G1A0
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
% Datum: 2020-08-06 19:02
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RRPRR12V1G1A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(4,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V1G1A0_inertia_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V1G1A0_inertia_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RRPRR12V1G1A0_inertia_para_pf_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR12V1G1A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRPRR12V1G1A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RRPRR12V1G1A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V1G1A0_inertia_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V1G1A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:01:31
% EndTime: 2020-08-06 19:01:33
% DurationCPUTime: 1.22s
% Computational Cost: add. (2409->223), mult. (3366->374), div. (330->6), fcn. (2355->18), ass. (0->183)
t575 = 2 * rSges(3,3);
t578 = (t575 + qJ(3,3)) * qJ(3,3);
t577 = (t575 + qJ(3,2)) * qJ(3,2);
t576 = (t575 + qJ(3,1)) * qJ(3,1);
t574 = m(2) * rSges(2,3);
t573 = (rSges(3,2) * m(3));
t494 = pkin(1) + rSges(3,1);
t572 = m(3) * t494;
t483 = sin(qJ(1,3));
t571 = t483 * pkin(4);
t485 = sin(qJ(1,2));
t570 = t485 * pkin(4);
t487 = sin(qJ(1,1));
t569 = t487 * pkin(4);
t482 = sin(qJ(2,3));
t568 = rSges(3,2) * t482;
t484 = sin(qJ(2,2));
t567 = rSges(3,2) * t484;
t486 = sin(qJ(2,1));
t566 = rSges(3,2) * t486;
t488 = cos(qJ(2,3));
t495 = pkin(1) + pkin(2);
t544 = t495 * t488;
t556 = t482 * qJ(3,3);
t456 = t544 + t556;
t489 = cos(qJ(1,3));
t469 = t489 * pkin(4);
t431 = t456 * t483 + t469;
t434 = t456 * t489 - t571;
t479 = legFrame(3,3);
t463 = sin(t479);
t466 = cos(t479);
t413 = -t463 * t431 + t466 * t434;
t497 = 1 / qJ(3,3);
t565 = t413 * t497;
t490 = cos(qJ(2,2));
t543 = t495 * t490;
t554 = t484 * qJ(3,2);
t457 = t543 + t554;
t491 = cos(qJ(1,2));
t470 = t491 * pkin(4);
t432 = t457 * t485 + t470;
t435 = t457 * t491 - t570;
t480 = legFrame(2,3);
t464 = sin(t480);
t467 = cos(t480);
t414 = -t464 * t432 + t467 * t435;
t499 = 1 / qJ(3,2);
t564 = t414 * t499;
t492 = cos(qJ(2,1));
t542 = t495 * t492;
t552 = t486 * qJ(3,1);
t458 = t542 + t552;
t493 = cos(qJ(1,1));
t471 = t493 * pkin(4);
t433 = t458 * t487 + t471;
t436 = t458 * t493 - t569;
t481 = legFrame(1,3);
t465 = sin(t481);
t468 = cos(t481);
t415 = -t465 * t433 + t468 * t436;
t501 = 1 / qJ(3,1);
t563 = t415 * t501;
t416 = t431 * t466 + t463 * t434;
t562 = t416 * t497;
t417 = t432 * t467 + t464 * t435;
t561 = t417 * t499;
t418 = t433 * t468 + t465 * t436;
t560 = t418 * t501;
t453 = -t488 * qJ(3,3) + t495 * t482;
t428 = (-t482 * t494 + t453) * t497 * m(3);
t559 = t428 * t497;
t454 = -t490 * qJ(3,2) + t495 * t484;
t429 = (-t484 * t494 + t454) * t499 * m(3);
t558 = t429 * t499;
t455 = -t492 * qJ(3,1) + t495 * t486;
t430 = (-t486 * t494 + t455) * t501 * m(3);
t557 = t430 * t501;
t555 = t482 * t497;
t553 = t484 * t499;
t551 = t486 * t501;
t550 = t488 * t497;
t549 = t490 * t499;
t548 = t492 * t501;
t547 = t494 * t497;
t546 = t494 * t499;
t545 = t494 * t501;
t541 = -Icges(2,1) - Icges(3,1);
t502 = rSges(3,3) ^ 2;
t540 = rSges(3,2) ^ 2 + t502;
t539 = rSges(3,3) - t494;
t538 = rSges(3,3) + t494;
t537 = m(3) * t547;
t536 = m(3) * t546;
t535 = m(3) * t545;
t437 = -t463 * t483 + t466 * t489;
t444 = t483 * t556 + t469;
t445 = t489 * t556 - t571;
t404 = t437 * t544 - t463 * t444 + t445 * t466;
t534 = t404 * t550;
t438 = t463 * t489 + t466 * t483;
t405 = t438 * t544 + t444 * t466 + t463 * t445;
t533 = t405 * t550;
t439 = -t464 * t485 + t467 * t491;
t446 = t485 * t554 + t470;
t447 = t491 * t554 - t570;
t406 = t439 * t543 - t464 * t446 + t447 * t467;
t532 = t406 * t549;
t440 = t464 * t491 + t467 * t485;
t407 = t440 * t543 + t446 * t467 + t464 * t447;
t531 = t407 * t549;
t441 = -t465 * t487 + t468 * t493;
t448 = t487 * t552 + t471;
t449 = t493 * t552 - t569;
t408 = t441 * t542 - t465 * t448 + t449 * t468;
t530 = t408 * t548;
t442 = t465 * t493 + t468 * t487;
t409 = t442 * t542 + t448 * t468 + t465 * t449;
t529 = t409 * t548;
t512 = pkin(1) ^ 2 + t502 + (2 * pkin(1) + rSges(3,1)) * rSges(3,1);
t504 = rSges(2,2) ^ 2;
t506 = rSges(2,1) ^ 2;
t516 = Icges(3,2) + Icges(2,3) + (t504 + t506) * m(2);
t425 = (t512 + t578) * m(3) + t516;
t419 = (t425 * t482 - t453 * t572) * t497;
t528 = t419 * t550;
t426 = (t512 + t577) * m(3) + t516;
t420 = (t426 * t484 - t454 * t572) * t499;
t527 = t420 * t549;
t427 = (t512 + t576) * m(3) + t516;
t421 = (t427 * t486 - t455 * t572) * t501;
t526 = t421 * t548;
t443 = rSges(2,1) * t574 + rSges(3,2) * t572 - Icges(3,4) - Icges(2,5);
t476 = rSges(3,3) + qJ(3,3);
t508 = -rSges(2,2) * t574 + Icges(2,6) - Icges(3,6);
t422 = ((t476 * t573) + t508) * t488 - t482 * t443;
t525 = t422 * t550;
t477 = rSges(3,3) + qJ(3,2);
t423 = ((t477 * t573) + t508) * t490 - t484 * t443;
t524 = t423 * t549;
t478 = rSges(3,3) + qJ(3,1);
t424 = ((t478 * t573) + t508) * t492 - t486 * t443;
t523 = t424 * t548;
t522 = t425 * t550;
t521 = t426 * t549;
t520 = t427 * t548;
t519 = t488 * t547;
t518 = t490 * t546;
t517 = t492 * t545;
t515 = t555 * t573;
t514 = t553 * t573;
t513 = t551 * t573;
t511 = Icges(1,3) + (rSges(2,3) ^ 2 + t504) * m(2) + (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) - t541;
t510 = (-t504 + t506) * m(2) + Icges(2,2) + Icges(3,3) + t541;
t509 = -m(2) * rSges(2,1) * rSges(2,2) + Icges(2,4) - Icges(3,5);
t452 = 0.1e1 / t458;
t451 = 0.1e1 / t457;
t450 = 0.1e1 / t456;
t412 = (t455 * t573 + t424) * t551;
t411 = (t454 * t573 + t423) * t553;
t410 = (t453 * t573 + t422) * t555;
t403 = (t540 + t576) * m(3) + (0.2e1 * t486 * (t478 * t572 + t509) + (-(qJ(3,1) + t538) * (qJ(3,1) + t539) * m(3) + t510) * t492) * t492 + t511;
t402 = (t540 + t577) * m(3) + (0.2e1 * t484 * (t477 * t572 + t509) + (-(qJ(3,2) + t538) * (qJ(3,2) + t539) * m(3) + t510) * t490) * t490 + t511;
t401 = (t540 + t578) * m(3) + (0.2e1 * t482 * (t476 * t572 + t509) + (-(qJ(3,3) + t538) * (qJ(3,3) + t539) * m(3) + t510) * t488) * t488 + t511;
t400 = (t563 + (-t408 * t517 - t442 * t566) * t452) * m(3);
t399 = (t560 + (-t409 * t517 + t441 * t566) * t452) * m(3);
t398 = (t564 + (-t406 * t518 - t440 * t567) * t451) * m(3);
t397 = (t561 + (-t407 * t518 + t439 * t567) * t451) * m(3);
t396 = (t565 + (-t404 * t519 - t438 * t568) * t450) * m(3);
t395 = (t562 + (-t405 * t519 + t437 * t568) * t450) * m(3);
t394 = -t415 * t535 + (t408 * t520 - t424 * t442) * t452;
t393 = -t418 * t535 + (t409 * t520 + t424 * t441) * t452;
t392 = -t414 * t536 + (t406 * t521 - t423 * t440) * t451;
t391 = -t417 * t536 + (t407 * t521 + t423 * t439) * t451;
t390 = -t413 * t537 + (t404 * t522 - t422 * t438) * t450;
t389 = -t416 * t537 + (t405 * t522 + t422 * t437) * t450;
t388 = t415 * t513 + (-t403 * t442 + t408 * t523) * t452;
t387 = t418 * t513 + (t403 * t441 + t409 * t523) * t452;
t386 = t414 * t514 + (-t402 * t440 + t406 * t524) * t451;
t385 = t417 * t514 + (t402 * t439 + t407 * t524) * t451;
t384 = t413 * t515 + (-t401 * t438 + t404 * t525) * t450;
t383 = t416 * t515 + (t401 * t437 + t405 * t525) * t450;
t1 = [t396 * t565 + t398 * t564 + t400 * t563 + m(4) + (-t388 * t442 + t394 * t530) * t452 + (-t386 * t440 + t392 * t532) * t451 + (-t384 * t438 + t390 * t534) * t450, t396 * t562 + t398 * t561 + t400 * t560 + (t388 * t441 + t394 * t529) * t452 + (t386 * t439 + t392 * t531) * t451 + (t384 * t437 + t390 * t533) * t450, (t394 * t486 + t400 * t455) * t501 + (t392 * t484 + t398 * t454) * t499 + (t390 * t482 + t396 * t453) * t497; t395 * t565 + t397 * t564 + t399 * t563 + (-t387 * t442 + t393 * t530) * t452 + (-t385 * t440 + t391 * t532) * t451 + (-t383 * t438 + t389 * t534) * t450, t395 * t562 + t397 * t561 + t399 * t560 + m(4) + (t387 * t441 + t393 * t529) * t452 + (t385 * t439 + t391 * t531) * t451 + (t383 * t437 + t389 * t533) * t450, (t393 * t486 + t399 * t455) * t501 + (t391 * t484 + t397 * t454) * t499 + (t389 * t482 + t395 * t453) * t497; t413 * t559 + t414 * t558 + t415 * t557 + (t408 * t526 - t412 * t442) * t452 + (t406 * t527 - t411 * t440) * t451 + (t404 * t528 - t410 * t438) * t450, t416 * t559 + t417 * t558 + t418 * t557 + (t409 * t526 + t412 * t441) * t452 + (t407 * t527 + t411 * t439) * t451 + (t405 * t528 + t410 * t437) * t450, m(4) + (t421 * t486 + t430 * t455) * t501 + (t420 * t484 + t429 * t454) * t499 + (t419 * t482 + t428 * t453) * t497;];
MX  = t1;
