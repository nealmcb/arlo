import { createLocation, createMemoryHistory } from 'history'
import { match as routerMatch } from 'react-router-dom'

type MatchParameter<Params> = { [K in keyof Params]?: string }

const generateUrl = <Params extends MatchParameter<Params>>(
  path: string,
  params: Params
): string => {
  let tempPath = path

  for (const param in params) {
    /* istanbul ignore else */
    if (params.hasOwnProperty(param)) {
      const value = params[param]
      tempPath = tempPath.replace(`:${param}`, value as NonNullable<
        typeof value
      >)
    }
  }

  return tempPath
}

export const routerTestProps = <Params extends MatchParameter<Params> = {}>(
  path: string,
  params: Params
) => {
  const match: routerMatch<Params> = {
    isExact: false,
    path,
    url: generateUrl(path, params),
    params,
  }
  const history = createMemoryHistory()
  const location = createLocation(match.url)

  return { history, location, match }
}

/** Credit to https://stackoverflow.com/a/56452779 for solution to mocking React Router props */

export const regexpEscape = (s: string) => {
  /* eslint-disable no-useless-escape */
  return s.replace(/[-\/\\^$*+?.()|[\]{}]/g, '\\$&')
}

export default {
  routerTestProps,
}
