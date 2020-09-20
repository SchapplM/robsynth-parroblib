% Calculate vector of centrifugal and coriolis load on the joints for
% P3PRRRR1G1P3A0
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
% Datum: 2020-03-09 20:34
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3PRRRR1G1P3A0_coriolisvec_para_pf_slag_vp2(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR1G1P3A0_coriolisvec_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR1G1P3A0_coriolisvec_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR1G1P3A0_coriolisvec_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRRR1G1P3A0_coriolisvec_para_pf_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR1G1P3A0_coriolisvec_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRRR1G1P3A0_coriolisvec_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3PRRRR1G1P3A0_coriolisvec_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR1G1P3A0_coriolisvec_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR1G1P3A0_coriolisvec_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 20:34:06
% EndTime: 2020-03-09 20:34:07
% DurationCPUTime: 0.74s
% Computational Cost: add. (1056->130), mult. (2880->276), div. (1806->10), fcn. (3591->18), ass. (0->120)
t496 = legFrame(3,3);
t477 = sin(t496);
t480 = cos(t496);
t507 = cos(qJ(3,3));
t501 = sin(qJ(3,3));
t508 = cos(qJ(2,3));
t541 = t501 * t508;
t456 = -t477 * t541 - t480 * t507;
t460 = -t477 * t507 + t480 * t541;
t513 = xDP(2);
t514 = xDP(1);
t515 = 0.1e1 / pkin(2);
t502 = sin(qJ(2,3));
t484 = 0.1e1 / t502;
t489 = 0.1e1 / t507 ^ 2;
t545 = t484 * t489;
t444 = (t456 * t514 + t460 * t513) * t515 * t545;
t488 = 0.1e1 / t507;
t450 = (t477 * t514 - t480 * t513) * t515 * t488;
t536 = t507 * t508;
t542 = t501 * t507;
t551 = t450 * t501;
t422 = ((t444 * t536 - t502 * t551) * t488 * t444 + (-t444 * t502 * t542 + t450 * t508) * t489 * t450) * t484;
t441 = t444 ^ 2;
t447 = t450 ^ 2;
t554 = t447 * t488;
t435 = (-t441 * t507 - t554) * t484 * pkin(2);
t471 = -Ifges(3,5) * t501 - Ifges(3,6) * t507;
t499 = Ifges(3,1) - Ifges(3,2);
t530 = t501 * t554;
t474 = mrSges(3,1) * t501 + mrSges(3,2) * t507;
t548 = t474 * t502;
t487 = t507 ^ 2;
t561 = 0.2e1 * t487;
t565 = t488 * (-Ifges(3,3) * t530 - t422 * t471 - t435 * t548 + t441 * (Ifges(3,4) * t561 + t499 * t542 - Ifges(3,4)));
t497 = legFrame(2,3);
t478 = sin(t497);
t481 = cos(t497);
t509 = cos(qJ(3,2));
t503 = sin(qJ(3,2));
t510 = cos(qJ(2,2));
t539 = t503 * t510;
t457 = -t478 * t539 - t481 * t509;
t462 = -t478 * t509 + t481 * t539;
t504 = sin(qJ(2,2));
t485 = 0.1e1 / t504;
t492 = 0.1e1 / t509 ^ 2;
t544 = t485 * t492;
t445 = (t457 * t514 + t462 * t513) * t515 * t544;
t491 = 0.1e1 / t509;
t451 = (t478 * t514 - t481 * t513) * t515 * t491;
t535 = t509 * t510;
t540 = t503 * t509;
t550 = t451 * t503;
t420 = ((t445 * t535 - t504 * t550) * t491 * t445 + (-t445 * t504 * t540 + t451 * t510) * t492 * t451) * t485;
t442 = t445 ^ 2;
t448 = t451 ^ 2;
t553 = t448 * t491;
t436 = (-t442 * t509 - t553) * t485 * pkin(2);
t472 = -Ifges(3,5) * t503 - Ifges(3,6) * t509;
t529 = t503 * t553;
t475 = mrSges(3,1) * t503 + mrSges(3,2) * t509;
t547 = t475 * t504;
t490 = t509 ^ 2;
t560 = 0.2e1 * t490;
t564 = t491 * (-Ifges(3,3) * t529 - t420 * t472 - t436 * t547 + t442 * (Ifges(3,4) * t560 + t499 * t540 - Ifges(3,4)));
t498 = legFrame(1,3);
t479 = sin(t498);
t482 = cos(t498);
t511 = cos(qJ(3,1));
t505 = sin(qJ(3,1));
t512 = cos(qJ(2,1));
t537 = t505 * t512;
t458 = -t479 * t537 - t482 * t511;
t464 = -t479 * t511 + t482 * t537;
t506 = sin(qJ(2,1));
t486 = 0.1e1 / t506;
t495 = 0.1e1 / t511 ^ 2;
t543 = t486 * t495;
t446 = (t458 * t514 + t464 * t513) * t515 * t543;
t494 = 0.1e1 / t511;
t452 = (t479 * t514 - t482 * t513) * t515 * t494;
t534 = t511 * t512;
t538 = t505 * t511;
t549 = t452 * t505;
t421 = ((t446 * t534 - t506 * t549) * t494 * t446 + (-t446 * t506 * t538 + t452 * t512) * t495 * t452) * t486;
t443 = t446 ^ 2;
t449 = t452 ^ 2;
t552 = t449 * t494;
t437 = (-t443 * t511 - t552) * t486 * pkin(2);
t473 = -Ifges(3,5) * t505 - Ifges(3,6) * t511;
t528 = t505 * t552;
t476 = mrSges(3,1) * t505 + mrSges(3,2) * t511;
t546 = t476 * t506;
t493 = t511 ^ 2;
t559 = 0.2e1 * t493;
t563 = t494 * (-Ifges(3,3) * t528 - t421 * t473 - t437 * t546 + t443 * (Ifges(3,4) * t559 + t499 * t538 - Ifges(3,4)));
t562 = -0.2e1 * Ifges(3,4);
t558 = Ifges(3,5) / 0.2e1;
t557 = -Ifges(3,6) / 0.2e1;
t500 = mrSges(2,2) - mrSges(3,3);
t556 = t500 / 0.2e1;
t555 = -Ifges(3,1) - Ifges(2,3);
t521 = -mrSges(3,1) * t507 + mrSges(3,2) * t501;
t453 = -(mrSges(2,1) - t521) * t508 + t502 * t500;
t483 = -m(1) - m(2) - m(3);
t533 = t453 * t422 - t530 * t548 + (-mrSges(2,1) * t441 + t521 * (t441 + t447)) * t502 - 0.2e1 * t444 * t508 * (t444 * t556 + t474 * t450) + t483 * t435;
t520 = -mrSges(3,1) * t509 + mrSges(3,2) * t503;
t454 = -(mrSges(2,1) - t520) * t510 + t504 * t500;
t532 = t454 * t420 - t529 * t547 + (-mrSges(2,1) * t442 + t520 * (t442 + t448)) * t504 - 0.2e1 * t445 * t510 * (t445 * t556 + t475 * t451) + t483 * t436;
t519 = -mrSges(3,1) * t511 + mrSges(3,2) * t505;
t455 = -(mrSges(2,1) - t519) * t512 + t506 * t500;
t531 = t455 * t421 - t528 * t546 + (-mrSges(2,1) * t443 + t519 * (t443 + t449)) * t506 - 0.2e1 * t446 * t512 * (t446 * t556 + t476 * t452) + t483 * t437;
t524 = t485 * t491 * t532;
t523 = t486 * t494 * t531;
t522 = t533 * t488 * t484;
t518 = (0.2e1 * ((t444 * t499 * t501 + t450 * t558) * t507 + t551 * t557 + (t561 - 0.1e1) * t444 * Ifges(3,4)) * t450 + t453 * t435 + (t487 * t499 + t542 * t562 + t555) * t422 - t471 * t530) * t545;
t517 = (0.2e1 * ((t445 * t499 * t503 + t451 * t558) * t509 + t550 * t557 + (t560 - 0.1e1) * t445 * Ifges(3,4)) * t451 + t454 * t436 + (t490 * t499 + t540 * t562 + t555) * t420 - t472 * t529) * t544;
t516 = (0.2e1 * ((t446 * t499 * t505 + t452 * t558) * t511 + t549 * t557 + (t559 - 0.1e1) * t446 * Ifges(3,4)) * t452 + t455 * t437 + (t493 * t499 + t538 * t562 + t555) * t421 - t473 * t528) * t543;
t1 = [(t479 * t505 + t482 * t534) * t523 + (t478 * t503 + t481 * t535) * t524 + (t477 * t501 + t480 * t536) * t522 + (t456 * t518 + t457 * t517 + t458 * t516 - t477 * t565 - t478 * t564 - t479 * t563) * t515; (t479 * t534 - t482 * t505) * t523 + (t478 * t535 - t481 * t503) * t524 + (t477 * t536 - t480 * t501) * t522 + (t460 * t518 + t462 * t517 + t464 * t516 + t480 * t565 + t481 * t564 + t482 * t563) * t515; t531 + t532 + t533;];
taucX  = t1;
