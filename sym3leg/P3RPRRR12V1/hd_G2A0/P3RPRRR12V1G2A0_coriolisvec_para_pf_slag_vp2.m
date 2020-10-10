% Calculate vector of centrifugal and coriolis load on the joints for
% P3RPRRR12V1G2A0
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
% Datum: 2020-08-06 18:25
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RPRRR12V1G2A0_coriolisvec_para_pf_slag_vp2(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR12V1G2A0_coriolisvec_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRRR12V1G2A0_coriolisvec_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR12V1G2A0_coriolisvec_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RPRRR12V1G2A0_coriolisvec_para_pf_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR12V1G2A0_coriolisvec_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RPRRR12V1G2A0_coriolisvec_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RPRRR12V1G2A0_coriolisvec_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR12V1G2A0_coriolisvec_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR12V1G2A0_coriolisvec_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:24:36
% EndTime: 2020-08-06 18:24:38
% DurationCPUTime: 1.55s
% Computational Cost: add. (6516->222), mult. (9585->374), div. (2883->7), fcn. (8313->18), ass. (0->152)
t593 = (-pkin(5) - pkin(6));
t500 = pkin(1) - t593;
t514 = sin(qJ(1,3));
t520 = cos(qJ(1,3));
t473 = qJ(2,3) * t520 - t500 * t514;
t510 = legFrame(3,2);
t494 = sin(t510);
t497 = cos(t510);
t513 = sin(qJ(3,3));
t519 = cos(qJ(3,3));
t572 = t519 * qJ(2,3);
t588 = pkin(3) * t520;
t589 = pkin(3) * t519;
t443 = (-t473 * t497 + t494 * t589) * t513 + (t519 - 0.1e1) * (t519 + 0.1e1) * t497 * t588 + t494 * t572;
t505 = t519 ^ 2;
t446 = (t473 * t494 + t497 * t589) * t513 + (-t505 + 0.1e1) * t494 * t588 + t497 * t572;
t483 = t513 * pkin(3) + qJ(2,3);
t470 = t514 * t483 + t500 * t520;
t480 = 0.1e1 / t483;
t502 = 0.1e1 / t513;
t525 = xDP(3);
t526 = xDP(2);
t527 = xDP(1);
t602 = 0.2e1 * (t470 * t525 + (t443 * t527 + t446 * t526) * t502) * t480;
t516 = sin(qJ(1,2));
t522 = cos(qJ(1,2));
t474 = qJ(2,2) * t522 - t500 * t516;
t511 = legFrame(2,2);
t495 = sin(t511);
t498 = cos(t511);
t515 = sin(qJ(3,2));
t521 = cos(qJ(3,2));
t571 = t521 * qJ(2,2);
t586 = pkin(3) * t522;
t587 = pkin(3) * t521;
t444 = (-t474 * t498 + t495 * t587) * t515 + (t521 - 0.1e1) * (t521 + 0.1e1) * t498 * t586 + t495 * t571;
t506 = t521 ^ 2;
t447 = (t474 * t495 + t498 * t587) * t515 + (-t506 + 0.1e1) * t495 * t586 + t498 * t571;
t484 = t515 * pkin(3) + qJ(2,2);
t471 = t516 * t484 + t500 * t522;
t481 = 0.1e1 / t484;
t503 = 0.1e1 / t515;
t601 = 0.2e1 * (t471 * t525 + (t444 * t527 + t447 * t526) * t503) * t481;
t518 = sin(qJ(1,1));
t524 = cos(qJ(1,1));
t475 = qJ(2,1) * t524 - t500 * t518;
t512 = legFrame(1,2);
t496 = sin(t512);
t499 = cos(t512);
t517 = sin(qJ(3,1));
t523 = cos(qJ(3,1));
t570 = t523 * qJ(2,1);
t584 = pkin(3) * t524;
t585 = pkin(3) * t523;
t445 = (-t475 * t499 + t496 * t585) * t517 + (t523 - 0.1e1) * (t523 + 0.1e1) * t499 * t584 + t496 * t570;
t507 = t523 ^ 2;
t448 = (t475 * t496 + t499 * t585) * t517 + (-t507 + 0.1e1) * t496 * t584 + t499 * t570;
t485 = t517 * pkin(3) + qJ(2,1);
t472 = t518 * t485 + t500 * t524;
t482 = 0.1e1 / t485;
t504 = 0.1e1 / t517;
t600 = 0.2e1 * (t472 * t525 + (t445 * t527 + t448 * t526) * t504) * t482;
t452 = (t520 * t525 + (-t494 * t526 + t497 * t527) * t514) * t480;
t536 = 0.1e1 / pkin(3);
t464 = (-t494 * t527 - t497 * t526) * t536 * t502;
t530 = qJ(2,3) ^ 2;
t535 = pkin(3) ^ 2;
t534 = pkin(5) ^ 2;
t537 = pkin(1) ^ 2;
t596 = -2 * pkin(5);
t542 = (2 * t593 * pkin(1)) - t534 - t535 - t537 + ((t596 - pkin(6)) * pkin(6));
t560 = t464 * t589;
t579 = qJ(2,3) * t513;
t561 = -0.2e1 * t579;
t575 = t500 * t452;
t592 = pkin(3) * t464;
t422 = (((t519 * t575 - t592) * t513 - qJ(2,3) * t464) * t502 * t592 + t575 * t602 + (t500 * t560 + (pkin(3) * t561 + t505 * t535 - t530 + t542) * t452) * t452) * t480;
t425 = (t602 + 0.2e1 * t560 - t575) * t452 * t480;
t461 = t464 ^ 2;
t528 = pkin(1) + pkin(5);
t489 = t528 * mrSges(3,2) - Ifges(3,6);
t490 = t528 * mrSges(3,1) - Ifges(3,5);
t467 = -t489 * t513 + t490 * t519;
t582 = (mrSges(3,3) - mrSges(2,2));
t479 = m(2) * pkin(1) + t528 * m(3) + t582;
t529 = m(2) + m(3);
t486 = t529 * qJ(2,3) + mrSges(2,3);
t509 = Ifges(3,1) - Ifges(3,2);
t538 = -m(3) * t534 + (mrSges(3,3) * t596) - t529 * t537 - Ifges(2,1) - Ifges(3,2) - Ifges(1,3);
t539 = (mrSges(3,1) * qJ(2,3) - t509 * t513) * t519;
t548 = -mrSges(3,1) * t513 - mrSges(3,2) * t519;
t554 = t461 * t502 * t519;
t557 = mrSges(3,2) * t579;
t583 = (m(3) * pkin(5) + t582) * pkin(1);
t595 = -0.2e1 * mrSges(2,3);
t569 = (mrSges(3,1) * t561 + qJ(2,3) * t595 - t509 * t505 - t529 * t530 + t538) * t425 + t479 * t422 + t467 * t554 + 0.2e1 * (-(mrSges(3,2) * qJ(2,3) - Ifges(3,4) * t513) * t519 - t583) * t425 + (t489 * t519 + t490 * t513) * t461 + ((t486 - t548) * t602 + (-0.2e1 * t557 + (-0.4e1 * t505 + 0.2e1) * Ifges(3,4) + 0.2e1 * t539) * t464) * t452;
t599 = t514 * t569;
t453 = (t522 * t525 + (-t495 * t526 + t498 * t527) * t516) * t481;
t465 = (-t495 * t527 - t498 * t526) * t536 * t503;
t531 = qJ(2,2) ^ 2;
t559 = t465 * t587;
t580 = qJ(2,2) * t515;
t562 = -0.2e1 * t580;
t574 = t500 * t453;
t591 = pkin(3) * t465;
t423 = (((t521 * t574 - t591) * t515 - qJ(2,2) * t465) * t503 * t591 + t574 * t601 + (t500 * t559 + (pkin(3) * t562 + t506 * t535 - t531 + t542) * t453) * t453) * t481;
t426 = (t601 + 0.2e1 * t559 - t574) * t453 * t481;
t462 = t465 ^ 2;
t468 = -t489 * t515 + t490 * t521;
t487 = t529 * qJ(2,2) + mrSges(2,3);
t540 = (mrSges(3,1) * qJ(2,2) - t509 * t515) * t521;
t547 = -mrSges(3,1) * t515 - mrSges(3,2) * t521;
t553 = t462 * t503 * t521;
t556 = mrSges(3,2) * t580;
t568 = (mrSges(3,1) * t562 + qJ(2,2) * t595 - t509 * t506 - t529 * t531 + t538) * t426 + t479 * t423 + t468 * t553 + 0.2e1 * (-(mrSges(3,2) * qJ(2,2) - Ifges(3,4) * t515) * t521 - t583) * t426 + (t489 * t521 + t490 * t515) * t462 + ((t487 - t547) * t601 + (-0.2e1 * t556 + (-0.4e1 * t506 + 0.2e1) * Ifges(3,4) + 0.2e1 * t540) * t465) * t453;
t598 = t516 * t568;
t454 = (t524 * t525 + (-t496 * t526 + t499 * t527) * t518) * t482;
t466 = (-t496 * t527 - t499 * t526) * t536 * t504;
t532 = qJ(2,1) ^ 2;
t558 = t466 * t585;
t581 = qJ(2,1) * t517;
t563 = -0.2e1 * t581;
t573 = t500 * t454;
t590 = pkin(3) * t466;
t424 = (((t523 * t573 - t590) * t517 - qJ(2,1) * t466) * t504 * t590 + t573 * t600 + (t500 * t558 + (pkin(3) * t563 + t507 * t535 - t532 + t542) * t454) * t454) * t482;
t427 = (t600 + 0.2e1 * t558 - t573) * t454 * t482;
t463 = t466 ^ 2;
t469 = -t489 * t517 + t490 * t523;
t488 = t529 * qJ(2,1) + mrSges(2,3);
t541 = (mrSges(3,1) * qJ(2,1) - t509 * t517) * t523;
t546 = -mrSges(3,1) * t517 - mrSges(3,2) * t523;
t552 = t463 * t504 * t523;
t555 = mrSges(3,2) * t581;
t567 = (mrSges(3,1) * t563 + qJ(2,1) * t595 - t509 * t507 - t529 * t532 + t538) * t427 + t479 * t424 + t469 * t552 + 0.2e1 * (-(mrSges(3,2) * qJ(2,1) - Ifges(3,4) * t517) * t523 - t583) * t427 + (t489 * t523 + t490 * t517) * t463 + ((t488 - t546) * t600 + (-0.2e1 * t555 + (-0.4e1 * t507 + 0.2e1) * Ifges(3,4) + 0.2e1 * t541) * t466) * t454;
t597 = t518 * t567;
t594 = -0.2e1 * Ifges(3,4);
t449 = t452 ^ 2;
t476 = -t519 * mrSges(3,1) + t513 * mrSges(3,2);
t566 = -t529 * t422 + t479 * t425 + t476 * t554 - t486 * t449 + t548 * (t449 + t461);
t450 = t453 ^ 2;
t477 = -t521 * mrSges(3,1) + t515 * mrSges(3,2);
t565 = -t529 * t423 + t479 * t426 + t477 * t553 - t487 * t450 + t547 * (t450 + t462);
t451 = t454 ^ 2;
t478 = -t523 * mrSges(3,1) + t517 * mrSges(3,2);
t564 = -t529 * t424 + t479 * t427 + t478 * t552 - t488 * t451 + t546 * (t451 + t463);
t551 = t566 * t502;
t550 = t565 * t503;
t549 = t564 * t504;
t545 = (Ifges(3,3) * t554 - t476 * t422 - t467 * t425 + t449 * (t505 * t594 + Ifges(3,4) + t539 - t557)) * t502;
t544 = (Ifges(3,3) * t553 - t477 * t423 - t468 * t426 + t450 * (t506 * t594 + Ifges(3,4) + t540 - t556)) * t503;
t543 = (Ifges(3,3) * t552 - t478 * t424 - t469 * t427 + t451 * (t507 * t594 + Ifges(3,4) + t541 - t555)) * t504;
t1 = [(t445 * t549 + t499 * t597) * t482 + (t444 * t550 + t498 * t598) * t481 + (t443 * t551 + t497 * t599) * t480 + (t494 * t545 + t495 * t544 + t496 * t543) * t536; (t448 * t549 - t496 * t597) * t482 + (t447 * t550 - t495 * t598) * t481 + (t446 * t551 - t494 * t599) * t480 + (t497 * t545 + t498 * t544 + t499 * t543) * t536; (t564 * t472 + t567 * t524) * t482 + (t565 * t471 + t568 * t522) * t481 + (t566 * t470 + t569 * t520) * t480;];
taucX  = t1;
