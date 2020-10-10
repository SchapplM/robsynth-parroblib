% Calculate vector of centrifugal and coriolis load on the joints for
% P3RPRRR12V1G1A0
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
% taucX [3x1]
%   forces required to compensate Coriolis and centrifugal joint torques
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:21
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RPRRR12V1G1A0_coriolisvec_para_pf_slag_vp1(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR12V1G1A0_coriolisvec_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRRR12V1G1A0_coriolisvec_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR12V1G1A0_coriolisvec_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RPRRR12V1G1A0_coriolisvec_para_pf_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR12V1G1A0_coriolisvec_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RPRRR12V1G1A0_coriolisvec_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RPRRR12V1G1A0_coriolisvec_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR12V1G1A0_coriolisvec_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR12V1G1A0_coriolisvec_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:21:20
% EndTime: 2020-08-06 18:21:21
% DurationCPUTime: 1.38s
% Computational Cost: add. (4614->224), mult. (6024->371), div. (1308->11), fcn. (5232->18), ass. (0->156)
t574 = 2 * m(3);
t573 = -2 * pkin(1);
t572 = 2 * qJ(2,1);
t571 = 2 * qJ(2,2);
t570 = 2 * qJ(2,3);
t481 = sin(qJ(3,3));
t464 = 0.1e1 / t481;
t487 = cos(qJ(3,3));
t555 = t464 * t487;
t465 = 0.1e1 / t481 ^ 2;
t495 = xDP(3);
t549 = t495 ^ 2 / pkin(3) ^ 2;
t533 = t465 * t549;
t483 = sin(qJ(3,2));
t467 = 0.1e1 / t483;
t489 = cos(qJ(3,2));
t553 = t467 * t489;
t468 = 0.1e1 / t483 ^ 2;
t531 = t468 * t549;
t485 = sin(qJ(3,1));
t470 = 0.1e1 / t485;
t491 = cos(qJ(3,1));
t551 = t470 * t491;
t471 = 0.1e1 / t485 ^ 2;
t529 = t471 * t549;
t569 = -0.2e1 * t487;
t568 = -0.2e1 * t489;
t567 = -0.2e1 * t491;
t566 = m(3) * rSges(3,2);
t493 = (rSges(3,3) + pkin(5));
t565 = m(3) * (t487 * rSges(3,1) - t481 * rSges(3,2));
t564 = m(3) * (t489 * rSges(3,1) - t483 * rSges(3,2));
t563 = m(3) * (t491 * rSges(3,1) - t485 * rSges(3,2));
t562 = (-pkin(1) - t493) * m(3);
t478 = legFrame(3,3);
t452 = sin(t478);
t455 = cos(t478);
t482 = sin(qJ(1,3));
t488 = cos(qJ(1,3));
t424 = -t452 * t482 + t455 * t488;
t427 = t452 * t488 + t455 * t482;
t446 = t481 * pkin(3) + qJ(2,3);
t443 = 0.1e1 / t446;
t496 = xDP(2);
t497 = xDP(1);
t406 = (t424 * t497 + t427 * t496) * t443;
t561 = t406 * t443;
t479 = legFrame(2,3);
t453 = sin(t479);
t456 = cos(t479);
t484 = sin(qJ(1,2));
t490 = cos(qJ(1,2));
t425 = -t453 * t484 + t456 * t490;
t428 = t453 * t490 + t456 * t484;
t447 = t483 * pkin(3) + qJ(2,2);
t444 = 0.1e1 / t447;
t407 = (t425 * t497 + t428 * t496) * t444;
t560 = t407 * t444;
t480 = legFrame(1,3);
t454 = sin(t480);
t457 = cos(t480);
t486 = sin(qJ(1,1));
t492 = cos(qJ(1,1));
t426 = -t454 * t486 + t457 * t492;
t429 = t454 * t492 + t457 * t486;
t448 = t485 * pkin(3) + qJ(2,1);
t445 = 0.1e1 / t448;
t408 = (t426 * t497 + t429 * t496) * t445;
t559 = t408 * t445;
t462 = pkin(1) + pkin(5) + pkin(6);
t558 = t462 * t406;
t557 = t462 * t407;
t556 = t462 * t408;
t554 = t464 * t495;
t552 = t467 * t495;
t550 = t470 * t495;
t507 = 0.1e1 / pkin(3);
t548 = t495 * t507;
t534 = t487 * t554;
t418 = t482 * t446 + t462 * t488;
t421 = -t446 * t488 + t462 * t482;
t409 = t418 * t455 - t452 * t421;
t412 = t452 * t418 + t421 * t455;
t541 = (t409 * t497 + t412 * t496) * t443;
t391 = t534 + t541;
t473 = t487 ^ 2;
t499 = qJ(2,3) ^ 2;
t506 = pkin(3) ^ 2;
t516 = pkin(3) * t462 * t548;
t526 = qJ(2,3) * t464 * t548;
t537 = -(t462 ^ 2) - t506;
t382 = ((t391 * t462 - t516 * t555) * t481 + ((t473 * t506 - t499 + t537) * t481 + (t473 - 0.1e1) * t570 * pkin(3)) * t406) * t464 * t561 + (t391 * t558 - ((t487 * t558 + t554) * t481 + t526) * t465 * t495) * t443;
t385 = (t391 - t534 + t541 - t558) * t561;
t437 = -rSges(3,2) * t562 - Icges(3,6);
t438 = -rSges(3,1) * t562 - Icges(3,5);
t415 = -t481 * t437 + t438 * t487;
t431 = -t562 + (pkin(1) - rSges(2,2)) * m(2);
t503 = rSges(3,2) ^ 2;
t505 = rSges(3,1) ^ 2;
t435 = (-t503 + t505) * m(3) - Icges(3,1) + Icges(3,2);
t439 = qJ(2,3) * m(3) + m(2) * (rSges(2,3) + qJ(2,3));
t450 = rSges(3,1) * t566 - Icges(3,4);
t517 = -(rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) - Icges(2,1) - Icges(3,2) - Icges(1,3);
t520 = -rSges(3,1) * t481 - rSges(3,2) * t487;
t509 = pkin(1) ^ 2;
t521 = rSges(2,3) ^ 2 + t509 + (t573 + rSges(2,2)) * rSges(2,2);
t522 = -t505 - t509 + ((t573 - t493) * t493);
t525 = t533 * t555;
t547 = (t435 * t473 + t450 * t481 * t569 - ((rSges(2,3) * t570 + t499 + t521) * m(2)) + (t520 * t570 - t499 + t522) * m(3) + t517) * t385 + t431 * t382 + t415 * t525 + (t437 * t487 + t438 * t481) * t533 + (0.2e1 * t391 * t439 + (t435 * t569 + (-0.4e1 * t473 + 0.2e1) * t464 * t450) * t548 + ((-rSges(3,1) * t526 + rSges(3,2) * t391) * t487 + (rSges(3,1) * t391 + rSges(3,2) * t526) * t481) * t574) * t406;
t532 = t489 * t552;
t419 = t484 * t447 + t462 * t490;
t422 = -t447 * t490 + t462 * t484;
t410 = t419 * t456 - t453 * t422;
t413 = t453 * t419 + t422 * t456;
t540 = (t410 * t497 + t413 * t496) * t444;
t392 = t532 + t540;
t474 = t489 ^ 2;
t500 = qJ(2,2) ^ 2;
t527 = qJ(2,2) * t467 * t548;
t383 = ((t392 * t462 - t516 * t553) * t483 + ((t474 * t506 - t500 + t537) * t483 + (t474 - 0.1e1) * t571 * pkin(3)) * t407) * t467 * t560 + (t392 * t557 - ((t489 * t557 + t552) * t483 + t527) * t468 * t495) * t444;
t386 = (t392 - t532 + t540 - t557) * t560;
t416 = -t483 * t437 + t438 * t489;
t440 = qJ(2,2) * m(3) + m(2) * (rSges(2,3) + qJ(2,2));
t519 = -rSges(3,1) * t483 - rSges(3,2) * t489;
t524 = t531 * t553;
t546 = (t435 * t474 + t450 * t483 * t568 - ((rSges(2,3) * t571 + t500 + t521) * m(2)) + (t519 * t571 - t500 + t522) * m(3) + t517) * t386 + t431 * t383 + t416 * t524 + (t437 * t489 + t438 * t483) * t531 + (0.2e1 * t392 * t440 + (t435 * t568 + (-0.4e1 * t474 + 0.2e1) * t467 * t450) * t548 + ((-rSges(3,1) * t527 + rSges(3,2) * t392) * t489 + (rSges(3,1) * t392 + rSges(3,2) * t527) * t483) * t574) * t407;
t530 = t491 * t550;
t420 = t486 * t448 + t462 * t492;
t423 = -t448 * t492 + t462 * t486;
t411 = t420 * t457 - t454 * t423;
t414 = t454 * t420 + t423 * t457;
t539 = (t411 * t497 + t414 * t496) * t445;
t393 = t530 + t539;
t475 = t491 ^ 2;
t501 = qJ(2,1) ^ 2;
t528 = qJ(2,1) * t470 * t548;
t384 = ((t393 * t462 - t516 * t551) * t485 + ((t475 * t506 - t501 + t537) * t485 + (t475 - 0.1e1) * t572 * pkin(3)) * t408) * t470 * t559 + (t393 * t556 - ((t491 * t556 + t550) * t485 + t528) * t471 * t495) * t445;
t387 = (t393 - t530 + t539 - t556) * t559;
t417 = -t485 * t437 + t438 * t491;
t441 = qJ(2,1) * m(3) + m(2) * (rSges(2,3) + qJ(2,1));
t518 = -rSges(3,1) * t485 - rSges(3,2) * t491;
t523 = t529 * t551;
t545 = (t435 * t475 + t450 * t485 * t567 - ((rSges(2,3) * t572 + t501 + t521) * m(2)) + (t518 * t572 - t501 + t522) * m(3) + t517) * t387 + t431 * t384 + t417 * t523 + (t437 * t491 + t438 * t485) * t529 + (0.2e1 * t393 * t441 + (t435 * t567 + (-0.4e1 * t475 + 0.2e1) * t470 * t450) * t548 + ((-rSges(3,1) * t528 + rSges(3,2) * t393) * t491 + (rSges(3,1) * t393 + rSges(3,2) * t528) * t485) * t574) * t408;
t403 = t406 ^ 2;
t498 = -m(2) - m(3);
t544 = t498 * t382 + t431 * t385 - t525 * t565 - t439 * t403 + t520 * (t403 + t533) * m(3);
t404 = t407 ^ 2;
t543 = t498 * t383 + t431 * t386 - t524 * t564 - t440 * t404 + t519 * (t404 + t531) * m(3);
t405 = t408 ^ 2;
t542 = t498 * t384 + t431 * t387 - t523 * t563 - t441 * t405 + t518 * (t405 + t529) * m(3);
t535 = -t566 / 0.2e1;
t538 = rSges(3,1) * t535 + Icges(3,4) / 0.2e1;
t536 = m(3) * rSges(3,1) / 0.2e1;
t442 = -(t503 + t505) * m(3) - Icges(3,3);
t430 = (t505 / 0.2e1 - t503 / 0.2e1) * m(3) - Icges(3,1) / 0.2e1 + Icges(3,2) / 0.2e1;
t1 = [(t542 * t411 + t545 * t426) * t445 + (t543 * t410 + t546 * t425) * t444 + (t544 * t409 + t547 * t424) * t443; (t542 * t414 + t545 * t429) * t445 + (t543 * t413 + t546 * t428) * t444 + (t544 * t412 + t547 * t427) * t443; t542 * t551 + t543 * t553 + t544 * t555 + ((-t442 * t523 + t384 * t563 - t417 * t387 + 0.2e1 * t405 * (t450 * t475 + (qJ(2,1) * t536 + t430 * t485) * t491 + t485 * qJ(2,1) * t535 + t538)) * t470 + (-t442 * t524 + t383 * t564 - t416 * t386 + 0.2e1 * t404 * (t450 * t474 + (qJ(2,2) * t536 + t430 * t483) * t489 + t483 * qJ(2,2) * t535 + t538)) * t467 + (-t442 * t525 + t382 * t565 - t415 * t385 + 0.2e1 * t403 * (t450 * t473 + (qJ(2,3) * t536 + t430 * t481) * t487 + t481 * qJ(2,3) * t535 + t538)) * t464) * t507;];
taucX  = t1;
