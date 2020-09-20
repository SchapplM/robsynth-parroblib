% Calculate vector of centrifugal and coriolis load on the joints for
% P3RPRR1G2A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% taucX [3x1]
%   forces required to compensate Coriolis and centrifugal joint torques
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:25
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RPRR1G2A0_coriolisvec_para_pf_slag_vp2(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRR1G2A0_coriolisvec_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRR1G2A0_coriolisvec_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRR1G2A0_coriolisvec_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRR1G2A0_coriolisvec_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRR1G2A0_coriolisvec_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RPRR1G2A0_coriolisvec_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RPRR1G2A0_coriolisvec_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRR1G2A0_coriolisvec_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRR1G2A0_coriolisvec_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:24:55
% EndTime: 2020-03-09 21:24:57
% DurationCPUTime: 1.26s
% Computational Cost: add. (14520->187), mult. (9417->253), div. (1818->7), fcn. (6798->59), ass. (0->135)
t544 = 2 * pkin(2);
t520 = 0.1e1 / pkin(3);
t548 = t520 / 0.2e1;
t506 = legFrame(1,2);
t543 = qJ(1,1) + pkin(7);
t483 = t506 + t543;
t477 = qJ(3,1) + t483;
t484 = -t506 + t543;
t478 = qJ(3,1) + t484;
t452 = sin(t477) + sin(t478);
t455 = -cos(t478) + cos(t477);
t501 = pkin(7) + qJ(3,1);
t509 = sin(qJ(3,1));
t460 = 0.1e1 / (pkin(1) * sin(t501) + t509 * pkin(2));
t487 = cos(qJ(1,1) + t501);
t514 = xDP(3);
t516 = xDP(1);
t549 = t516 / 0.2e1;
t515 = xDP(2);
t550 = t515 / 0.2e1;
t419 = (t452 * t549 + t455 * t550 + t487 * t514) * t460;
t496 = qJ(1,1) + t506;
t497 = qJ(1,1) - t506;
t428 = -t455 * pkin(3) + (cos(t484) - cos(t483)) * pkin(2) + (cos(t497) - cos(t496)) * pkin(1);
t540 = t515 * t548;
t422 = t428 * t460 * t540;
t449 = -pkin(3) * t487 - pkin(2) * cos(t543) - cos(qJ(1,1)) * pkin(1);
t547 = t514 * t520;
t443 = t449 * t460 * t547;
t539 = t516 * t548;
t551 = (t452 * pkin(3) + (sin(t483) + sin(t484)) * pkin(2) + (sin(t496) + sin(t497)) * pkin(1)) * t460;
t530 = t539 * t551;
t407 = -t530 / 0.2e1 + t422 / 0.2e1 + t443 / 0.2e1 + t419;
t413 = t422 + t443 - t530;
t410 = t413 + t419;
t503 = cos(pkin(7));
t488 = t503 * pkin(1) + pkin(2);
t491 = cos(t501);
t521 = pkin(2) ^ 2;
t522 = pkin(1) ^ 2;
t498 = t521 + t522;
t512 = cos(qJ(3,1));
t519 = pkin(3) ^ 2;
t545 = 0.2e1 * pkin(1);
t546 = pkin(3) * t544;
t554 = t410 * t413;
t558 = pkin(2) * t503;
t502 = sin(pkin(7));
t559 = pkin(1) * t502;
t401 = (t407 * t512 * t546 + t410 * t519 + t419 * t498 + (pkin(3) * t407 * t491 + t419 * t558) * t545) * t520 * t460 * t419 + (t512 * t488 - t509 * t559 + pkin(3)) / (t509 * t488 + t512 * t559) * t554;
t404 = (-pkin(3) * t554 + (-t410 * pkin(3) + (-pkin(1) * t491 - t512 * pkin(2)) * t419) * t419) * t460;
t456 = mrSges(3,1) * t559 + t488 * mrSges(3,2);
t457 = t488 * mrSges(3,1) - mrSges(3,2) * t559;
t431 = t456 * t509 - t457 * t512 - Ifges(3,3);
t434 = t456 * t512 + t509 * t457;
t523 = -(m(3) * t521) - (m(2) + m(3)) * t522 - Ifges(1,3) - Ifges(2,3) - Ifges(3,3);
t533 = -t512 * mrSges(3,1) + mrSges(3,2) * t509;
t557 = m(3) * pkin(2) + mrSges(2,1);
t524 = (-0.2e1 * t407 * t413 * t434 + (t533 * t544 + (-(-t533 + t557) * t503 + (mrSges(3,1) * t509 + t512 * mrSges(3,2) + mrSges(2,2)) * t502) * t545 + t523) * t404 + t431 * t401) * t460;
t562 = t524 / 0.2e1;
t505 = legFrame(2,2);
t542 = qJ(1,2) + pkin(7);
t481 = t505 + t542;
t475 = qJ(3,2) + t481;
t482 = -t505 + t542;
t476 = qJ(3,2) + t482;
t451 = sin(t475) + sin(t476);
t454 = -cos(t476) + cos(t475);
t500 = pkin(7) + qJ(3,2);
t508 = sin(qJ(3,2));
t459 = 0.1e1 / (pkin(1) * sin(t500) + t508 * pkin(2));
t486 = cos(qJ(1,2) + t500);
t418 = (t451 * t549 + t454 * t550 + t486 * t514) * t459;
t494 = qJ(1,2) + t505;
t495 = qJ(1,2) - t505;
t427 = -t454 * pkin(3) + (cos(t482) - cos(t481)) * pkin(2) + (cos(t495) - cos(t494)) * pkin(1);
t421 = t427 * t459 * t540;
t448 = -pkin(3) * t486 - pkin(2) * cos(t542) - cos(qJ(1,2)) * pkin(1);
t442 = t448 * t459 * t547;
t552 = (t451 * pkin(3) + (sin(t481) + sin(t482)) * pkin(2) + (sin(t494) + sin(t495)) * pkin(1)) * t459;
t531 = t539 * t552;
t406 = -t531 / 0.2e1 + t421 / 0.2e1 + t442 / 0.2e1 + t418;
t412 = t421 + t442 - t531;
t409 = t412 + t418;
t490 = cos(t500);
t511 = cos(qJ(3,2));
t555 = t409 * t412;
t400 = (t406 * t511 * t546 + t409 * t519 + t418 * t498 + (pkin(3) * t406 * t490 + t418 * t558) * t545) * t520 * t459 * t418 + (t488 * t511 - t508 * t559 + pkin(3)) / (t488 * t508 + t511 * t559) * t555;
t403 = (-pkin(3) * t555 + (-pkin(3) * t409 + (-pkin(1) * t490 - t511 * pkin(2)) * t418) * t418) * t459;
t430 = t456 * t508 - t457 * t511 - Ifges(3,3);
t433 = t456 * t511 + t508 * t457;
t534 = -t511 * mrSges(3,1) + mrSges(3,2) * t508;
t525 = (-0.2e1 * t406 * t412 * t433 + (t534 * t544 + (-(-t534 + t557) * t503 + (mrSges(3,1) * t508 + t511 * mrSges(3,2) + mrSges(2,2)) * t502) * t545 + t523) * t403 + t430 * t400) * t459;
t561 = t525 / 0.2e1;
t504 = legFrame(3,2);
t541 = qJ(1,3) + pkin(7);
t479 = t504 + t541;
t473 = qJ(3,3) + t479;
t480 = -t504 + t541;
t474 = qJ(3,3) + t480;
t450 = sin(t473) + sin(t474);
t453 = -cos(t474) + cos(t473);
t499 = pkin(7) + qJ(3,3);
t507 = sin(qJ(3,3));
t458 = 0.1e1 / (pkin(1) * sin(t499) + t507 * pkin(2));
t485 = cos(qJ(1,3) + t499);
t417 = (t450 * t549 + t453 * t550 + t485 * t514) * t458;
t492 = qJ(1,3) + t504;
t493 = qJ(1,3) - t504;
t426 = -t453 * pkin(3) + (cos(t480) - cos(t479)) * pkin(2) + (cos(t493) - cos(t492)) * pkin(1);
t420 = t426 * t458 * t540;
t447 = -pkin(3) * t485 - pkin(2) * cos(t541) - cos(qJ(1,3)) * pkin(1);
t441 = t447 * t458 * t547;
t553 = (t450 * pkin(3) + (sin(t479) + sin(t480)) * pkin(2) + (sin(t492) + sin(t493)) * pkin(1)) * t458;
t532 = t539 * t553;
t405 = -t532 / 0.2e1 + t420 / 0.2e1 + t441 / 0.2e1 + t417;
t411 = t420 + t441 - t532;
t408 = t411 + t417;
t489 = cos(t499);
t510 = cos(qJ(3,3));
t556 = t408 * t411;
t399 = (t405 * t510 * t546 + t408 * t519 + t417 * t498 + (pkin(3) * t405 * t489 + t417 * t558) * t545) * t520 * t458 * t417 + (t488 * t510 - t507 * t559 + pkin(3)) / (t488 * t507 + t510 * t559) * t556;
t402 = (-pkin(3) * t556 + (-pkin(3) * t408 + (-pkin(1) * t489 - t510 * pkin(2)) * t417) * t417) * t458;
t429 = t456 * t507 - t457 * t510 - Ifges(3,3);
t432 = t456 * t510 + t507 * t457;
t535 = -t510 * mrSges(3,1) + mrSges(3,2) * t507;
t526 = (-0.2e1 * t405 * t411 * t432 + (t535 * t544 + (-(-t535 + t557) * t503 + (mrSges(3,1) * t507 + t510 * mrSges(3,2) + mrSges(2,2)) * t502) * t545 + t523) * t402 + t429 * t399) * t458;
t560 = t526 / 0.2e1;
t538 = t417 ^ 2 * t432 - Ifges(3,3) * t399 + t429 * t402;
t537 = t418 ^ 2 * t433 - Ifges(3,3) * t400 + t430 * t403;
t536 = t419 ^ 2 * t434 - Ifges(3,3) * t401 + t431 * t404;
t529 = t538 * t458;
t528 = t537 * t459;
t527 = t536 * t460;
t1 = [t452 * t562 + t451 * t561 + t450 * t560 + (-t536 * t551 - t537 * t552 - t538 * t553) * t548; t455 * t562 + t454 * t561 + t453 * t560 + (t426 * t529 + t427 * t528 + t428 * t527) * t548; t487 * t524 + t486 * t525 + t485 * t526 + (t447 * t529 + t448 * t528 + t449 * t527) * t520;];
taucX  = t1;
