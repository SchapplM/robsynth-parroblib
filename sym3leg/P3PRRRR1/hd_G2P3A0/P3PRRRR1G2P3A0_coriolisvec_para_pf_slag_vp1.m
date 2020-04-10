% Calculate vector of centrifugal and coriolis load on the joints for
% P3PRRRR1G2P3A0
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
% taucX [3x1]
%   forces required to compensate Coriolis and centrifugal joint torques
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:16
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3PRRRR1G2P3A0_coriolisvec_para_pf_slag_vp1(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR1G2P3A0_coriolisvec_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR1G2P3A0_coriolisvec_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR1G2P3A0_coriolisvec_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRRR1G2P3A0_coriolisvec_para_pf_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR1G2P3A0_coriolisvec_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRRR1G2P3A0_coriolisvec_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3PRRRR1G2P3A0_coriolisvec_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR1G2P3A0_coriolisvec_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR1G2P3A0_coriolisvec_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:16:24
% EndTime: 2020-03-09 21:16:25
% DurationCPUTime: 0.94s
% Computational Cost: add. (1224->136), mult. (3525->283), div. (2598->10), fcn. (3717->24), ass. (0->125)
t485 = sin(qJ(2,1));
t467 = 0.1e1 / t485;
t490 = cos(qJ(3,1));
t475 = 0.1e1 / t490;
t495 = xDP(3);
t503 = 0.1e1 / pkin(2);
t479 = legFrame(1,2);
t460 = sin(t479);
t463 = cos(t479);
t496 = xDP(2);
t497 = xDP(1);
t504 = t460 * t496 - t463 * t497;
t476 = 0.1e1 / t490 ^ 2;
t484 = sin(qJ(3,1));
t491 = cos(qJ(2,1));
t526 = t476 * t484 * t491;
t425 = (-t475 * t495 - t504 * t526) * t503 * t467;
t437 = t504 * t503 * t475;
t540 = t485 * t490;
t546 = t425 * t491;
t410 = ((-t437 * t485 * t484 + t490 * t546) * t475 * t425 + (-t484 * t425 * t540 + t491 * t437) * t476 * t437) * t467;
t422 = t425 ^ 2;
t434 = t437 ^ 2;
t543 = t434 * t475;
t416 = (-t422 * t490 - t543) * t467 * pkin(2);
t552 = rSges(3,3) * m(3);
t455 = rSges(3,2) * t552 - Icges(3,6);
t456 = rSges(3,1) * t552 - Icges(3,5);
t440 = t455 * t490 + t456 * t484;
t501 = rSges(3,2) ^ 2;
t502 = rSges(3,1) ^ 2;
t539 = (t501 + t502);
t453 = -t539 * m(3) - Icges(3,3);
t457 = m(3) * rSges(3,1) * rSges(3,2) - Icges(3,4);
t474 = t490 ^ 2;
t454 = m(2) * rSges(2,2) - t552;
t494 = m(2) * rSges(2,1);
t559 = (rSges(3,1) * t490 - rSges(3,2) * t484) * m(3);
t431 = -(t494 + t559) * t491 + t485 * t454;
t500 = 0.2e1 * qJ(3,1);
t510 = -(rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) - (2 * rSges(3,3) ^ 2 + t539) * m(3) / 0.2e1 - Icges(3,2) / 0.2e1 - Icges(3,1) / 0.2e1 - Icges(2,3);
t451 = (-t501 + t502) * m(3) - Icges(3,1) + Icges(3,2);
t529 = t451 * t484 * t490;
t532 = t484 * t543;
t549 = t456 / 0.4e1;
t550 = -t455 / 0.4e1;
t551 = -t451 / 0.2e1;
t507 = t467 * (0.4e1 * ((t484 * t550 + t490 * t549) * t437 + (t529 / 0.2e1 + (t474 - 0.1e1 / 0.2e1) * t457) * t425) * t437 - t431 * t416 - (cos(t500) * t551 + t457 * sin(t500) + t510) * t410 + t440 * t532);
t449 = t484 * rSges(3,1) + t490 * rSges(3,2);
t535 = m(3) * t449 * t485;
t553 = 0.2e1 * t457;
t562 = t507 * t526 + (t440 * t410 + t416 * t535 + t422 * (t474 * t553 - t457 + t529) - t453 * t532) * t475;
t483 = sin(qJ(2,2));
t466 = 0.1e1 / t483;
t488 = cos(qJ(3,2));
t472 = 0.1e1 / t488;
t478 = legFrame(2,2);
t459 = sin(t478);
t462 = cos(t478);
t505 = t459 * t496 - t462 * t497;
t473 = 0.1e1 / t488 ^ 2;
t482 = sin(qJ(3,2));
t489 = cos(qJ(2,2));
t527 = t473 * t482 * t489;
t424 = (-t472 * t495 - t505 * t527) * t503 * t466;
t436 = t505 * t503 * t472;
t541 = t483 * t488;
t547 = t424 * t489;
t409 = ((-t436 * t483 * t482 + t488 * t547) * t472 * t424 + (-t482 * t424 * t541 + t489 * t436) * t473 * t436) * t466;
t421 = t424 ^ 2;
t433 = t436 ^ 2;
t544 = t433 * t472;
t415 = (-t421 * t488 - t544) * t466 * pkin(2);
t439 = t455 * t488 + t456 * t482;
t471 = t488 ^ 2;
t558 = (rSges(3,1) * t488 - rSges(3,2) * t482) * m(3);
t430 = -(t494 + t558) * t489 + t483 * t454;
t499 = 0.2e1 * qJ(3,2);
t530 = t451 * t482 * t488;
t533 = t482 * t544;
t508 = t466 * (0.4e1 * ((t482 * t550 + t488 * t549) * t436 + (t530 / 0.2e1 + (t471 - 0.1e1 / 0.2e1) * t457) * t424) * t436 - t430 * t415 - (cos(t499) * t551 + t457 * sin(t499) + t510) * t409 + t439 * t533);
t448 = t482 * rSges(3,1) + t488 * rSges(3,2);
t536 = m(3) * t448 * t483;
t561 = t508 * t527 + (t439 * t409 + t415 * t536 + t421 * (t471 * t553 - t457 + t530) - t453 * t533) * t472;
t481 = sin(qJ(2,3));
t465 = 0.1e1 / t481;
t486 = cos(qJ(3,3));
t469 = 0.1e1 / t486;
t477 = legFrame(3,2);
t458 = sin(t477);
t461 = cos(t477);
t506 = t458 * t496 - t461 * t497;
t470 = 0.1e1 / t486 ^ 2;
t480 = sin(qJ(3,3));
t487 = cos(qJ(2,3));
t528 = t470 * t480 * t487;
t423 = (-t469 * t495 - t506 * t528) * t503 * t465;
t435 = t506 * t503 * t469;
t542 = t481 * t486;
t548 = t423 * t487;
t408 = ((-t435 * t481 * t480 + t486 * t548) * t469 * t423 + (-t480 * t423 * t542 + t487 * t435) * t470 * t435) * t465;
t420 = t423 ^ 2;
t432 = t435 ^ 2;
t545 = t432 * t469;
t414 = (-t420 * t486 - t545) * t465 * pkin(2);
t438 = t455 * t486 + t456 * t480;
t468 = t486 ^ 2;
t557 = (rSges(3,1) * t486 - rSges(3,2) * t480) * m(3);
t429 = -(t494 + t557) * t487 + t481 * t454;
t498 = 0.2e1 * qJ(3,3);
t531 = t451 * t480 * t486;
t534 = t480 * t545;
t509 = t465 * (0.4e1 * ((t480 * t550 + t486 * t549) * t435 + (t531 / 0.2e1 + (t468 - 0.1e1 / 0.2e1) * t457) * t423) * t435 - t429 * t414 - (cos(t498) * t551 + t457 * sin(t498) + t510) * t408 + t438 * t534);
t447 = t480 * rSges(3,1) + t486 * rSges(3,2);
t537 = m(3) * t447 * t481;
t560 = t509 * t528 + (t438 * t408 + t414 * t537 + t420 * (t468 * t553 - t457 + t531) - t453 * t534) * t469;
t538 = 0.2e1 * m(3);
t464 = -m(1) - m(2) - m(3);
t525 = t465 * (t429 * t408 + t464 * t414 - t534 * t537 + (-t420 * t494 - (t420 + t432) * t557) * t481 - (t447 * t435 * t538 + t423 * t454) * t548);
t524 = t466 * (t430 * t409 + t464 * t415 - t533 * t536 + (-t421 * t494 - (t421 + t433) * t558) * t483 - (t448 * t436 * t538 + t424 * t454) * t547);
t523 = t467 * (t431 * t410 + t464 * t416 - t532 * t535 + (-t422 * t494 - (t422 + t434) * t559) * t485 - (t449 * t437 * t538 + t425 * t454) * t546);
t516 = t469 * t525;
t515 = t472 * t524;
t514 = t475 * t523;
t1 = [(t460 * t540 - t463 * t484) * t514 + (t459 * t541 - t462 * t482) * t515 + (t458 * t542 - t461 * t480) * t516 + (-t560 * t461 - t561 * t462 - t562 * t463) * t503; (t460 * t484 + t463 * t540) * t514 + (t459 * t482 + t462 * t541) * t515 + (t458 * t480 + t461 * t542) * t516 + (t560 * t458 + t561 * t459 + t562 * t460) * t503; t491 * t523 + t489 * t524 + t487 * t525 + (t469 * t509 + t472 * t508 + t475 * t507) * t503;];
taucX  = t1;
