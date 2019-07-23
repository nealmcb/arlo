import React from 'react'
import styled from 'styled-components'
import { toast } from 'react-toastify'
import FormSection, {
  FormSectionLabel,
  FormSectionDescription,
} from '../Form/FormSection'
import FormWrapper from '../Form/FormWrapper'
import FormButton from '../Form/FormButton'
import FormField from '../Form/FormField'
import FormButtonBar from '../Form/FormButtonBar'
import { api } from '../utilities'

const InputSection = styled.div`
  display: block;
  margin-top: 25px;
  width: 100%;
  font-size: 0.4em;
`

const InputLabel = styled.label`
  display: inline-block;
`

const InlineInput = styled.div`
  display: flex;
  flex-direction: row;
  justify-content: space-between;
  margin-bottom: 10px;
  width: 50%;
`

interface Props {
  audit: any
  isLoading: boolean
  setIsLoading: (isLoading: boolean) => void
  updateAudit: () => void
}

const CalculateRiskMeasurmeent = (props: Props) => {
  const { audit, isLoading, setIsLoading, updateAudit } = props
  if (!audit) {
    return <></>
  }

  const downloadBallotRetrievalList = (id: number, e: any) => {
    e.preventDefault()
    const jurisdictionID: string = audit.jurisdictions[0].id
    window.open(`/jurisdiction/${jurisdictionID}/${id}/retrieval-list`)
  }

  const downloadAuditReport = async (i: number, round: any, evt: any) => {
    evt.preventDefault()
    try {
      window.open(`/audit/report`)
      updateAudit()
    } catch (err) {
      toast.error(err.message)
    }
  }

  const calculateRiskMeasurement = async (data: any, evt: any) => {
    evt.preventDefault()
    const { id, candidateOne, candidateTwo } = data
    try {
      const jurisdictionID: string = audit.jurisdictions[0].id
      const body: any = {
        contests: [
          {
            id: 'contest-1',
            results: {
              'candidate-1': Number(candidateOne),
              'candidate-2': Number(candidateTwo),
            },
          },
        ],
      }

      setIsLoading(true)
      await api(`/jurisdiction/${jurisdictionID}/${id}/results`, {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify(body),
      })
      updateAudit()
    } catch (err) {
      toast.error(err.message)
    }
  }

  return audit.rounds.map((v: any, i: number) => {
    const round: number = i + 1
    const contest: any = v.contests.length > 0 ? v.contests[0] : undefined
    let candidateOne = ''
    let candidateTwo = ''
    const showCalculateButton =
      i + 1 === audit.rounds.length &&
      contest &&
      contest.endMeasurements &&
      !contest.endMeasurements.isComplete
    /* eslint-disable react/no-array-index-key */
    return (
      <React.Fragment key={i}>
        <FormWrapper title={`Round ${i + 1}`}>
          <FormSection
            label={`Ballot Retrieval List \n
              ${contest ? `${contest.sampleSize} Ballots` : ''}`}
          >
            {/*<SectionLabel>
              Ballot Retrieval List \n
              {contest ? `${contest.sampleSize} Ballots` : ''}
            </SectionLabel>*/}
            <FormButton
              onClick={(e: React.MouseEvent) =>
                downloadBallotRetrievalList(round, e)
              }
              inline
            >
              Download Ballot Retrieval List for Round {i + 1}
            </FormButton>
            <FormSectionLabel>Audited Results: Round {round}</FormSectionLabel>
            <FormSectionDescription>
              Enter the number of votes recorded for each candidate/choice in
              the audited ballots for Round {i + 1}
            </FormSectionDescription>
            <form>
              <InputSection>
                <InlineInput
                  onChange={(e: any) => (candidateOne = e.target.value)}
                >
                  <InputLabel>{audit.contests[0].choices[0].name}</InputLabel>
                  <FormField />
                </InlineInput>
                <InlineInput
                  onChange={(e: any) => (candidateTwo = e.target.value)}
                >
                  <InputLabel>{audit.contests[0].choices[1].name}</InputLabel>
                  <FormField />
                </InlineInput>
              </InputSection>
            </form>
          </FormSection>
          {isLoading && <p>Loading...</p>}
          {showCalculateButton && !isLoading && (
            <FormButtonBar>
              <FormButton
                onClick={(e: any) =>
                  calculateRiskMeasurement(
                    {
                      id: round,
                      round: v,
                      candidateOne,
                      candidateTwo,
                      roundIndex: i,
                    },
                    e
                  )
                }
              >
                Calculate Risk Measurement
              </FormButton>
            </FormButtonBar>
          )}
          {contest &&
            contest.endMeasurements.pvalue &&
            contest.endMeasurements.isComplete && (
              <FormSection>
                <FormSectionLabel>
                  Audit Status:{' '}
                  {contest.endMeasurements.isComplete
                    ? 'COMPLETE'
                    : 'INCOMPLETE'}
                </FormSectionLabel>
                <InputSection>
                  <InlineInput>
                    <InputLabel>Risk Limit: </InputLabel>
                    {audit.riskLimit}%
                  </InlineInput>
                  <InlineInput>
                    <InputLabel>P-value: </InputLabel>{' '}
                    {contest.endMeasurements.pvalue}
                  </InlineInput>
                </InputSection>
                {/* {Form 3} */}
                {contest.endMeasurements.isComplete && (
                  <FormButton
                    onClick={(e: React.MouseEvent) =>
                      downloadAuditReport(i, v, e)
                    }
                    size="sm"
                    inline
                  >
                    Download Audit Report
                  </FormButton>
                )}
              </FormSection>
            )}
        </FormWrapper>
      </React.Fragment>
    )
  })
}

export default React.memo(CalculateRiskMeasurmeent)